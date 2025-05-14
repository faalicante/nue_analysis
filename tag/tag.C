#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TText.h"
#include "TNtuple.h"
#include "TList.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include <iostream>
#include <filesystem>
#include <map>
#include <vector>
#include <stdexcept>
#include <cmath>

void printMemoryInfo() {
  ProcInfo_t procInfo;  // Declare a ProcInfo_t structure
  gSystem->GetProcInfo(&procInfo);  // Pass its pointer to GetProcInfo

  Long64_t mem = procInfo.fMemResident;  // Access the memory resident field
  std::cout << "Memory used: " << mem/1024 << " MB" << std::endl;
}

// Parameters
const int binSize    = 50;   //um
const int radius     = 300;  //um
const int ntag = 5;
const int bkg = 400;
int xMin, xMax, yMin, yMax, xCell, yCell;
int nPlates;
TString path;
TString opath;

void setPath(int data, TString* path, TString* opath, int cell, int* xMin, int* xMax, int* yMin, int* yMax, int* nPlates, int* xCell, int* yCell) {
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
      *path = "/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b000121/cells";
      *opath = "/eos/experiment/sndlhc/users/falicant/RUN1/b121/tag_trk";
      *nPlates = 57;
      *xCell = cell%18 +1;
      *yCell = cell/18 +1;
      *xMin = *xCell*10000 - 5000;
      *xMax = *xCell*10000 + 5000;
      *yMin = *yCell*10000 - 5000;
      *yMax = *yCell*10000 + 5000;
    }
}

void set_range(TH3F &h3, TH2F &h2) {
  std::cout << xMin << " " << xMax << std::endl;
  h3.GetXaxis()->SetRangeUser(xMin,xMax);
  h3.GetYaxis()->SetRangeUser(yMin,yMax);
  h2.GetXaxis()->SetRangeUser(xMin,xMax);
  h2.GetYaxis()->SetRangeUser(yMin,yMax);
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
    for(int iix = ix-r0; iix<=ix+r0; iix++) {//to fix remove the circle not square
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

void getEntriesInEllipse(TH2F &h2, TEllipse &el, int *entries, const int bkg) {
  int nBinsX = h2.GetNbinsX();
  int nBinsY = h2.GetNbinsY();
  float x0 = el.GetX1();
  float y0 = el.GetY1();
  float r1 = el.GetR1();
  *entries = 0;
  int bkg_plate = (int)round((double)bkg/nPlates);

  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
      float x = h2.GetXaxis()->GetBinCenter(i);
      float y = h2.GetYaxis()->GetBinCenter(j);
      float distance = (pow(x-x0,2)+pow(y-y0,2))/pow(r1,2);
      if (distance <= 1) {
        int binContent = h2.GetBinContent(i, j);
        if (binContent >= bkg_plate) {
          *entries = *entries + binContent - bkg_plate;
        }
      }
    }
  }  
  if (*entries<0) *entries = 0;
}

void count_bins(TH2F &h2, TObjArray &peaks, int plate, TH1F **h_long, const int bkg) {
  int np = peaks.GetEntries();
  int entries;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    getEntriesInEllipse(h2, *el, &entries, bkg);
    h_long[i]->SetBinContent(plate, entries);
  }
}

void makePlots(int cell, TCanvas *c, int np, int npmax, TH1F *h_long, int *maxPeak, int *maxPlate) {
  int idx = (np%3) +1;
  if (idx==1) c->Clear("D");
  c->cd(idx)->SetGrid(1,0);
  *maxPeak = h_long->GetMaximum(),
  *maxPlate = h_long->GetMaximumBin();
  h_long->SetLineColor(1);
  h_long->SetLineWidth(2);
  h_long->Draw("hist");
  if(np == 2)                 c->Print(Form("%s/longitudinal_xz_%i.pdf(", opath.Data(), cell), "pdf");
  else if(np == npmax-1)      c->Print(Form("%s/longitudinal_xz_%i.pdf)", opath.Data(), cell), "pdf");
  else if(idx == 3 && np >2 ) c->Print(Form("%s/longitudinal_xz_%i.pdf", opath.Data(), cell), "pdf");
}

void findStart(TH1F* h_long, int *firstPlate, int *lastPlate) {
  float entries = h_long->Integral();
  *firstPlate = h_long->FindFirstBinAbove();
  *lastPlate = h_long->FindLastBinAbove();
  // std::cout << entries << " " << entries * 0.01 << " " << *firstPlate << std::endl;
}

void makeNtuple(int cell, TH1F **h_long, TObjArray &peaks, int *ranks) {
  TString outputFileName = TString::Format("%s/peaks_%i.root", opath.Data(), cell);
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  TNtuple *ntuple = new TNtuple("showers","tagged showers","cell:tag:x:y:start:end:peak:maxplate:nseg:rankbin");
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 1500);
  c2->Divide(1,3);
  int np = peaks.GetEntries();
  int entries, firstPlate, lastPlate, maxPeak, maxPlate;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    float x = el->GetX1();
    float y = el->GetY1();
    findStart(h_long[i], &firstPlate, &lastPlate);
    makePlots(cell, c2, i, np, h_long[i], &maxPeak, &maxPlate);
    int nseg = h_long[i]->Integral(firstPlate, lastPlate);
    ntuple->Fill(cell, i+1, x, y, firstPlate, lastPlate, maxPeak, maxPlate, nseg, ranks[i]);
  }

  outputFile->Write();
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <data>" << argv[1] << " <cell>" << argv[2] << std::endl;
      // data = {0: muon, 1: nue, 2: data}
      // cell for MC is brick
      return 1;
  }
  int data = std::atoi(argv[1]);
  int cell = std::atoi(argv[2]);
  
  TStopwatch stopWatch;
  stopWatch.Start();
  
  gStyle->SetOptStat(0);
  TH2::AddDirectory(false);
  
  setPath(data, &path, &opath, cell, &xMin, &xMax, &yMin, &yMax, &nPlates, &xCell, &yCell);
  
  std::cout << "Cell " << cell << " x " << xCell << " y " << yCell << std::endl;

  TString file =  TString::Format("%s/cell_%i0_%i0/b000021/b000021.0.0.0.trk.root", path.Data(), xCell, yCell);
  TFile* f = TFile::Open(file);
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file: " << file << std::endl;
    return 1;
  }

  TH3F *h3 = (TH3F *)(f->Get("XYPseg"));
  TH2F *h2 = (TH2F *)(f->Get("XYseg"));
  set_range(*h3,*h2);

  h2->Smooth();
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  c->SetGrid();
  gStyle->SetOptStat(0);
  h2->Draw("colz");
  c->Update();
  c->Print(Form("%s/sh_%i.gif+180", opath.Data(), cell));
  TObjArray peaks;
  TObjArray txt;
  int ranks[ntag];
  get_peaks(*h2,peaks,txt,ntag,ranks,bkg);
  // h2->GetZaxis()->SetRangeUser(bkg, h2->GetMaximum());
  c->Update();
  c->Print(Form("%s/sh_%i.gif+180", opath.Data(), cell));
  drawEllipse(peaks,txt, kBlack);
  c->Update();
  c->Print(Form("%s/sh_%i.gif+180", opath.Data(), cell));

  
  TH1F *h_long[ntag];
  for(int i=0; i<ntag; i++) {
    h_long[i] = new TH1F(Form("h_long_%i", i+1),Form("Cluster %i;plate;%%segments", i+1), 60, 1, 61);
  }
  for(int p=1; p<=nPlates; p++) { 
    // printf("Plate %i\n", p);
    h3->GetZaxis()->SetRange(p+1, p+1);
    TH2F *h = (TH2F*)(h3->Project3D("yx"));
    h->Smooth();
    h->SetTitle(Form("Plate %i", p));
    h->Draw("colz");
    drawEllipse(peaks,txt, kBlack);
    count_bins(*h, peaks, p, &h_long[0], bkg);
    // h->GetZaxis()->SetRangeUser((int)round((double)bkg/nPlates), h->GetMaximum());
    c->Update();
    c->Print(Form("%s/sh_%i.gif+12", opath.Data(), cell));
  }
  c->Print(Form("%s/sh_%i.gif++", opath.Data(), cell));

  makeNtuple(cell, &h_long[0], peaks,ranks);

  std::cout << "Time: " << stopWatch.RealTime() << std::endl;
  printMemoryInfo();

  return 0;
}
