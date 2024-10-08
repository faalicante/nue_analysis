#include "TFile.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TText.h"
#include "TNtuple.h"
#include "TStyle.h"
#include <iostream>

const char *path = "/Users/fabioali/cernbox";
// const char *path = "/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b000022/shift";

int getMax(TH2F &h2, TObjArray &peaks, TObjArray &txt) {
  Int_t MaxBin = h2.GetMaximumBin();
  int rankbin = h2.GetMaximum();
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
  int r0=10;
  for(int iix = ix-r0; iix<=ix+r0; iix++) //to fix remove the circle not square
    for(int iiy = iy-r0; iiy<=iy+r0; iiy++)
      h2.SetBinContent(iix,iiy,0);
  return rankbin;
}

void get_peaks(TH2F &h2, TObjArray &peaks, TObjArray &txt, int npmax, int ranks[npmax]) {
  TH2F *h2new = (TH2F*)h2.Clone("get_peaks");
  for(int i=0; i<npmax; i++){
    int rankbin = getMax(*h2new, peaks, txt);
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

void getEntriesInEllipse(TH2F &h2, TEllipse &el, int *signif_bins, float *entries, float *bkg) {
  int nBinsX = h2.GetNbinsX();
  int nBinsY = h2.GetNbinsY();
  float x0 = el.GetX1();
  float y0 = el.GetY1();
  float r1 = el.GetR1();
  *signif_bins = 0;
  *entries = 0;

  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
      float x = h2.GetXaxis()->GetBinCenter(i);
      float y = h2.GetYaxis()->GetBinCenter(j);
      float distance = (pow(x-x0,2)+pow(y-y0,2))/pow(r1,2);
      if (distance <= 1) {
        int binContent = h2.GetBinContent(i, j);
        *entries = *entries + binContent - *bkg;
        }
      }
    }
  // }
  
  if (*entries<0) *entries = 0;
}

void count_bins(int npmax, TH2F &h2, TObjArray &peaks, int plate, TH1F *h_long[npmax]) {
  int np = peaks.GetEntries();
  int signif_bins;
  float bkg = 0;
  float entries;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    // eval_bkg(h2, *el, &bkg);
    getEntriesInEllipse(h2, *el, &signif_bins, &entries, &bkg);
    if (entries<0) entries = 0;
    h_long[i]->SetBinContent(plate, entries);
  }
}

void makePlots(int combination, TCanvas *c, int np, int npmax, TH1F* h_long, int *maxPeak, int *maxPlate) {
  int idx = (np%3) +1;
  if (idx==1) c->Clear("D");
  c->cd(idx)->SetGrid(1,0);
  *maxPeak = h_long->GetMaximum(),
  *maxPlate = h_long->GetMaximumBin();
  h_long->SetLineColor(1);
  h_long->SetLineWidth(2);
  h_long->Draw("hist");
  if(np == 2)                 c->Print(Form("%s/longitudinal_xz_%i.pdf(", path, combination), "pdf");
  else if(np == npmax-1)      c->Print(Form("%s/longitudinal_xz_%i.pdf)", path, combination), "pdf");
  else if(idx == 3 && np >2 ) c->Print(Form("%s/longitudinal_xz_%i.pdf", path, combination), "pdf");
}

void findStart(TH1F* h_long, int *firstPlate, int *lastPlate) { //add second axis?
  float entries = h_long->Integral();
  *firstPlate = h_long->FindFirstBinAbove(0.01 * entries);
  *lastPlate = h_long->FindLastBinAbove(0.01 * entries);
}

void makeNtuple(int combination, int np, TH1F *h_long[np], TObjArray &peaks, int ranks[np]) {
  TString outputFileName = TString::Format("%s/shift_full/peaks_shift_%d.root", path, combination);
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  TFile *output = new TFile(Form("%s/peaks_shift_%i.root", path,combination),"RECREATE");
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

  output->Write();
  output->Close();
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <combination>" << std::endl;
      return 1;
  }
  int combination = std::atoi(argv[1]);
  TString file =  TString::Format("%s/shift_full/histo_shifts_%d.root", path, combination);
  TFile* f = TFile::Open(file);

  int ntag = 100;
  TH2F *h2 = (TH2F *)(f->Get("XYseg"));
  // h2->Smooth();
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  c->SetGrid();
  gStyle->SetOptStat(0);
  h2->Draw("colz");
  c->Update();
  c->Print(Form("%s/sh_%i.gif+180", path, combination));
  
  TObjArray peaks;
  TObjArray txt;
  int ranks[ntag];
  get_peaks(*h2,peaks,txt,ntag,ranks);
  drawEllipse(peaks,txt, kBlack);
  c->Update();
  c->Print(Form("%s/sh_%i.gif+180", path, combination));

  
  TH1F *h_long[ntag];
  for(int i=0; i<ntag; i++) {
    h_long[i] = new TH1F(Form("h_long_%i", i+1),Form("Cluster %i;plate;%%segments", i+1), 60, 1, 61);
  }

  const int nplates = 60; //change nplates for mc
  for(int p=1; p<=nplates; p++) { 
    printf("Plate %i\n", p);
    TH2F *h = (TH2F *)(f->Get(Form("XYPseg_%i",p)));
    // h->Smooth();
    h->SetTitle(Form("Plate %i",p));
    h->Draw("colz");

    drawEllipse(peaks,txt, kBlack);
    count_bins(ntag, *h, peaks, p, &h_long[0]);
    c->Update();
    c->Print(Form("%s/sh_%i.gif+12", path, combination));
  }
  c->Print(Form("%s/sh_%i.gif++", path, combination));

  makeNtuple(combination, ntag, &h_long[0], peaks,ranks);

  return 0;
}
