#include "TFile.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TText.h"
#include "TNtuple.h"
#include "TStyle.h"
#include <iostream>

int nPlates;
TString path;

void setPath(int data, TString* path, int cell, int* nPlates) {
    if (data==0) {
        // *path = TString::Format("/Users/fabioali/cernbox/shift_muon/%i", cell);
        *path = TString::Format("/eos/experiment/sndlhc/users/falicant/shift_muon/%i", cell);
        *nPlates = 60;
    }
    else if (data==1) {
        *path = TString::Format("/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b%05i", cell);
        *nPlates = 60;
    }
    else if (data==2) {
        // *path = TString("/Users/fabioali/cernbox");
        *path = TString::Format("/eos/experiment/sndlhc/users/falicant/RUN1/b121/shift/%i", cell);
        *nPlates = 57;
    }
}

int getMax(TH2F &h2, TObjArray &peaks, TObjArray &txt, const int bkg) {
  int rankbin = h2.GetMaximum();
  if (rankbin>bkg) {
    Int_t MaxBin = h2.GetMaximumBin();
    Int_t ix,iy,iz;
    int radius = 300;
    h2.GetBinXYZ(MaxBin, ix, iy, iz);
    float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
    float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
    TEllipse  *el = new TEllipse(x,y,radius,radius);
    el->SetFillStyle(0);
    peaks.Add(el);
    TText  *t = new TText(x,y+300,Form("%d",peaks.GetEntries()));
    t->SetTextSize(0.03);
    txt.Add(t);
    int r0=6;
    for(int iix = ix-r0; iix<=ix+r0; iix++) {
      for(int iiy = iy-r0; iiy<=iy+r0; iiy++) {
        float distance = (pow(iix-ix,2)+pow(iiy-iy,2))/pow(r0,2);
        if (distance <= 1) h2.SetBinContent(iix,iiy,0);
      }
    }  
  return rankbin;
  }
  return 0;
}

void get_peaks(TH2F &h2, TObjArray &peaks, TObjArray &txt, int npmax, int *ranks, const int bkg) {
  TH2F *h2new = (TH2F*)h2.Clone("get_peaks");
  for(int i=0; i<npmax; i++){
    int rankbin = getMax(*h2new, peaks, txt, bkg);
    ranks[i] = rankbin;
  }
}

void drawEllipse(TObjArray &peaks, TObjArray &txt, int col) {
  int np = peaks.GetEntries();
  TText t(0,0,"a");
  for(int j=0; j<np; j++) {
    TEllipse *l = ((TEllipse*)(peaks.At(j)));
    l->SetLineColor(col);
    l->Draw();
    if(col!=kWhite) ((TText*)(txt.At(j)))->Draw();
  }
}

void getEntriesInEllipse(TH2F &h2, TEllipse &el, float *entries) {
  int nBinsX = h2.GetNbinsX();
  int nBinsY = h2.GetNbinsY();
  float x0 = el.GetX1();
  float y0 = el.GetY1();
  float r1 = el.GetR1();
  *entries = 0;

  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
      float x = h2.GetXaxis()->GetBinCenter(i);
      float y = h2.GetYaxis()->GetBinCenter(j);
      float distance = (pow(x-x0,2)+pow(y-y0,2))/pow(r1,2);
      if (distance <= 1) {
        int binContent = h2.GetBinContent(i, j);
        *if(binContent>8) *entries = *entries + binContent;
        }
      }
    }  
  if (*entries<0) *entries = 0;
}

void count_bins(TH2F &h2, TObjArray &peaks, int plate, TH1F **h_long, int data) {
  int np = peaks.GetEntries();
  float entries;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    getEntriesInEllipse(h2, *el, &entries);
    if (entries<0) entries = 0;
    h_long[i]->SetBinContent(plate, entries);
  }
}

void makePlots(int combination, TCanvas *c, int np, int npmax, TH1F *h_long, int *maxPeak, int *maxPlate) {
  int idx = (np%3) +1;
  if (idx==1) c->Clear("D");
  c->cd(idx)->SetGrid(1,0);
  *maxPeak = h_long->GetMaximum(),
  *maxPlate = h_long->GetMaximumBin();
  h_long->SetLineColor(1);
  h_long->SetLineWidth(2);
  h_long->Draw("hist");
  if(np == 2)                 c->Print(Form("%s/longitudinal_xz_%i.pdf(", path.Data(), combination), "pdf");
  else if(np == npmax-1)      c->Print(Form("%s/longitudinal_xz_%i.pdf)", path.Data(), combination), "pdf");
  else if(idx == 3 && np >2 ) c->Print(Form("%s/longitudinal_xz_%i.pdf", path.Data(), combination), "pdf");
}

void findStart(TH1F* h_long, int *firstPlate, int *lastPlate) { //add second axis?
  float entries = h_long->Integral();
  *firstPlate = h_long->FindFirstBinAbove(800);
  *lastPlate = h_long->FindLastBinAbove(800);
}

void makeNtuple(int combination, int np, TH1F **h_long, TObjArray &peaks, int *ranks) {
  TString outputFileName = TString::Format("%s/peaks_shift_%d.root", path.Data(), combination);
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  TNtuple *ntuple = new TNtuple("showers","tagged showers","combination:tag:x:y:start:end:peak:maxplate:nseg:rankbin");
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 1500);
  c2->Divide(1,3);
  
  float entries;
  float width;
  int firstPlate, lastPlate, maxPeak, maxPlate;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    float x = el->GetX1();
    float y = el->GetY1();
    findStart(h_long[i], &firstPlate, &lastPlate);
    makePlots(combination, c2, i, np, h_long[i], &maxPeak, &maxPlate);
    int nseg = h_long[i]->Integral(firstPlate, lastPlate);
    ntuple->Fill(combination, i+1, x, y, firstPlate, lastPlate, maxPeak, maxPlate, nseg, ranks[i]);
  }

  outputFile->Write();
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
      std::cerr << "Usage: " << argv[0] << " <data>" << argv[1] << " <combination>" << argv[2] << " <cell>" << argv[3] << std::endl;
      // data = {0: muon, 1: nue, 2: data}
      // cell for MC is brick
      return 1;
  }
  int data = std::atoi(argv[1]);
  int combination = std::atoi(argv[2]);
  int cell = std::atoi(argv[3]);
  std::cout << "Combination " << combination << std::endl;
  
  setPath(data, &path, cell, &nPlates);

  TString file =  TString::Format("%s/histo_shifts_%i.root", path.Data(), combination);
  TFile* f = TFile::Open(file);

  int ntag = 10;
  const int bkg = 500;
  TH2F *h2 = (TH2F *)(f->Get("XYseg"));
  h2->Smooth();
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  // c->SetGrid();
  // gStyle->SetOptStat(0);
  h2->Draw("colz");
  // c->Update();
  // c->Print(Form("%s/sh_%i.gif+180", path.Data(), combination));
  TObjArray peaks;
  TObjArray txt;
  int ranks[ntag];
  get_peaks(*h2,peaks,txt,ntag,ranks);
  h2->GetZaxis()->SetRangeUser(bkg, ranks[0]);
  // c->Update();
  // c->Print(Form("%s/sh_%i.gif+180", path.Data(), combination));
  drawEllipse(peaks,txt, kBlack);
  // c->Update();
  // c->Print(Form("%s/sh_%i.gif+180", path.Data(), combination));

  
  TH1F *h_long[ntag];
  for(int i=0; i<ntag; i++) {
    h_long[i] = new TH1F(Form("h_long_%i", i+1),Form("Cluster %i;plate;%%segments", i+1), 60, 1, 61);
  }

  for(int p=1; p<=nPlates; p++) { 
    // printf("Plate %i\n", p);
    TH2F *h = (TH2F *)(f->Get(Form("XYPseg_%i",p)));
    h->Smooth();
    h->SetTitle(Form("Plate %i",p));
    h->Draw("colz");

    drawEllipse(peaks,txt, kBlack);
    count_bins(ntag, *h, peaks, p, &h_long[0], data);
    // c->Update();
    // c->Print(Form("%s/sh_%i.gif+12", path.Data(), combination));
  }
  // c->Print(Form("%s/sh_%i.gif++", path.Data(), combination));

  makeNtuple(combination, ntag, &h_long[0], peaks,ranks);

  return 0;
}
