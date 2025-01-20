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
const int shiftRange = 50;   //mrad
const int shiftStep  = 2;    //mrad
const int stepZ      = 1350; //um
const int radius     = 300;  //um
const int ntag = 20;
const int bkg = 340;
int xMin, xMax, yMin, yMax, xBin, yBin, xLow, yLow;
int nPlates;
TString path;
TString opath;

int getBrick(int brick) {
    int bricks1[] = {11, 21, 31, 41, 51};
    int bricks2[] = {12, 22, 32, 42, 52};
    int bricks3[] = {13, 23, 33, 43, 53};
    int bricks4[] = {14, 24, 34, 44, 54};
    for (int i = 0; i < 5; i++) {
        if (brick == bricks1[i]) return 1;
        else if (brick == bricks2[i]) return 2;
        else if (brick == bricks3[i]) return 3;
        else if (brick == bricks4[i]) return 4;
    }
    return 0;
}

void setRange(int data, TString* path, TString* opath, int cell, int* xMin, int* xMax, int* yMin, int* yMax, int* xBin, int* yBin, int* nPlates, int* xLow, int* yLow) {
    if (data==0) {
        // *path = "/Users/fabioali/cernbox";
        *path = TString::Format("/eos/experiment/sndlhc/users/falicant/simulations/muon1E5_simsndlhc/b%06i", cell);
        *opath = *path;
        *nPlates = 60;
        *xMin = 284000;
        *xMax = 304000;
        *yMin = 79000;
        *yMax = 99000;
    }
    else if (data==1) {
        *path = TString::Format("/eos/experiment/sndlhc/users/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b%06i", cell);
        *opath = *path;
        *nPlates = 60;
        int group = getBrick(cell);
        switch (group) {
            case 1:
                *xMin = 200000;
                *yMin = 0;
                break;
            case 2:
                *xMin = 0;
                *yMin = 0;
                break;
            case 3:
                *xMin = 200000;
                *yMin = 200000;
                break;
            case 4:
                *xMin = 0;
                *yMin = 200000;
                break;
        }
        *xMax = *xMin + 200000;
        *yMax = *yMin + 200000;
    }
    else if (data==2) {
        // *path = "/Users/fabioali/cernbox/test_shift/trk";
        *path = "/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b000121/cells";
        *opath = TString::Format("/eos/experiment/sndlhc/users/falicant/RUN1/b121/shift/trk/%i", cell);
        *nPlates = 57;
        *xLow = cell % 18;
        *yLow = cell / 18;
        *xMin = *xLow*10000 + 4000;
        *xMax = *xLow*10000 + 16000;
        *yMin = *yLow*10000 + 4000;
        *yMax = *yLow*10000 + 16000;
    }
    *xBin = int((*xMax - *xMin) / binSize);
    *yBin = int((*yMax - *yMin) / binSize);
}

void openFiles(int cell, TFile* f[9]) {
    int idx = 0;
    for (int yCell = yLow-1; yCell <= yLow+1; yCell ++) {
        for (int xCell = xLow-1; xCell <= xLow+1; xCell++) {
            if (xCell < 1 || xCell > 18 || yCell < 1 || yCell > 18) f[idx] = nullptr;
            else {
                TString histFile = TString::Format("%s/cell_%i0_%i0_1x1cm/b000021/b000021.0.0.0.trk.root", path.Data(), xCell+1, yCell+1);
                f[idx] = TFile::Open(histFile);
            }
            idx++;   
        }
    }
}

TH2F* projectHist(TFile* f, int plate) {
    TString histName = TString::Format("XYseg_%d", plate);
    TH3F* h3 = (TH3F*)(f->Get("XYPseg"));
    h3->GetZaxis()->SetRange(plate+1,plate+1);
    TH2F* h2 = (TH2F*)(h3->Project3D("yx"));
    return h2;
}

TH2F* matrixCells(TFile* f[9], int plate, double shiftX, double shiftY) {
    TString histName = TString::Format("XYseg_%d", plate);
    TH2F* hm = new TH2F(histName, histName, xBin, xMin, xMax, yBin, yMin, yMax);
    TH2F* h2;
    for (int i = 0; i < 9; i++) {
        if (f[i] == nullptr) continue;
        h2 = projectHist(f[i], plate);
        for (int xBin = 1; xBin <= h2->GetNbinsX(); ++xBin) {
            double xCenter = h2->GetXaxis()->GetBinCenter(xBin) + shiftX;
            if (xCenter <= xMax && xCenter >= xMin) {
                for (int yBin = 1; yBin <= h2->GetNbinsY(); ++yBin) {
                    double yCenter = h2->GetYaxis()->GetBinCenter(yBin) + shiftY;
                    if (yCenter <= yMax && yCenter >= yMin) {
                        double content = h2->GetBinContent(xBin, yBin);
                        if (content) {
                            int xBinNew = hm->GetXaxis()->FindBin(xCenter);
                            int yBinNew = hm->GetYaxis()->FindBin(yCenter);
                            if (content > hm->GetBinContent(xBinNew, yBinNew)) hm->SetBinContent(xBinNew, yBinNew, content);
                        }
                    }
                }
            }
        }
        delete h2;
    }
    return hm;
}

TH2F* stackHist(int data, int combination, int cell, TH2F **hm) {
    TFile* f, *ff[9];

    if (data==0) {
        TString fileName = TString::Format("%s/hist_XY_muon.root", path.Data());
        f = TFile::Open(fileName);
        // h3 = (TH3F*)(f->Get("XYPseg"));
    }
    else if (data==1) {
        TString fileName = TString::Format("%s/hist_XYP_nue.root", path.Data());
        f = TFile::Open(fileName);
    }
    else if (data==2) {
        openFiles(cell, &ff[0]);
    }

    double shiftTX = (combination / (shiftRange+1)) * shiftStep - shiftRange;
    double shiftTY = (combination % (shiftRange+1)) * shiftStep - shiftRange;
    // std::cout << combination << " " << shiftRange+1 << " " << shiftStep << " " << shiftTX << " " << shiftTY << std::endl;
    // TH2F* h2;
    TH2F* hComb = new TH2F("XYseg", "XYseg", xBin, xMin, xMax, yBin, yMin, yMax);
    for (int layer = 0; layer < nPlates; ++layer) {
        int plate = layer + 1;
        std::cout << "Plate " << plate << std::endl;
        double shiftX = shiftTX / 1000.0 * stepZ * layer;
        double shiftY = shiftTY / 1000.0 * stepZ * layer;

        TString histName = TString::Format("XYseg_%d", plate);
        if (data ==0 || data==1) {
            TH2F *h2 = (TH2F*)(f->Get(histName.Data()));
            hm[layer] = matrixCells(&ff[0], plate, shiftX, shiftY);

        }
        else if (data==2) {
            // h2 = matrixCells(&ff[0], plate);

            hm[layer] = matrixCells(&ff[0], plate, shiftX, shiftY);
        }
        // delete h2;
        hComb->Add(hm[layer]);
    }
    if (data!=2) f->Close();
    else {
        for (int i = 0; i < 9; i++) {
            if (ff[i] != nullptr) ff[i]->Close();
        }
    }
    return hComb;
}

int getMax(TH2F &h2, TObjArray &peaks, TObjArray &txt, const int bkg) {
  int rankbin = h2.GetMaximum();
  if (rankbin>bkg) {
    Int_t MaxBin = h2.GetMaximumBin();
    Int_t ix,iy,iz;
    h2.GetBinXYZ(MaxBin, ix, iy, iz);
    float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
    float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
    TEllipse  *el = new TEllipse(x,y,radius,radius);
    el->SetFillStyle(0);
    peaks.Add(el);
    TText  *t = new TText(x,y+300,Form("%d",peaks.GetEntries()));
    t->SetTextSize(0.03);
    txt.Add(t);
    int r0=(int)round((double)radius/binSize);
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

void getEntriesInEllipse(TH2F &h2, TEllipse &el, int *entries, const int bkg) {
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
                if(binContent>= static_cast<int>(std::ceil((double)bkg/nPlates))) {
                    *entries = *entries + binContent;
                }
            }
        } 
    }
}

void count_bins(TH2F *h2, TObjArray &peaks, int plate, TH1F **h_long, const int bkg) {
  int np = peaks.GetEntries();
  int entries;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    getEntriesInEllipse(*h2, *el, &entries, bkg);
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
  if(np == 2)                 c->Print(Form("%s/longitudinal_xz_%i.pdf(", opath.Data(), combination), "pdf");
  else if(np == npmax-1)      c->Print(Form("%s/longitudinal_xz_%i.pdf)", opath.Data(), combination), "pdf");
  else if(idx == 3 && np >2 ) c->Print(Form("%s/longitudinal_xz_%i.pdf", opath.Data(), combination), "pdf");
}

void findStart(TH1F* h_long, int *firstPlate, int *lastPlate) {
  float entries = h_long->Integral();
  *firstPlate = h_long->FindFirstBinAbove(entries * 0.01);
  *lastPlate = h_long->FindLastBinAbove(entries * 0.01);
//   std::cout << entries << " " << entries * 0.01 << " " << *firstPlate << std::endl;
}

void makeNtuple(int combination, int cell, TH1F **h_long, TObjArray &peaks, int *ranks) {
  TString outputFileName = TString::Format("%s/peaks_%i.root", opath.Data(), combination);
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  TNtuple *ntuple = new TNtuple("showers","tagged showers","cell:combination:tag:x:y:start:end:peak:maxplate:nseg:rankbin");
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 1500);
  c2->Divide(1,3);
  int np = peaks.GetEntries();
  int entries, firstPlate, lastPlate, maxPeak, maxPlate;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    float x = el->GetX1();
    float y = el->GetY1();
    findStart(h_long[i], &firstPlate, &lastPlate);
    makePlots(combination, c2, i, np, h_long[i], &maxPeak, &maxPlate);
    int nseg = h_long[i]->Integral(firstPlate, lastPlate);
    ntuple->Fill(cell, combination, i+1, x, y, firstPlate, lastPlate, maxPeak, maxPlate, nseg, ranks[i]);
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

    TStopwatch stopWatch;
    stopWatch.Start();

    setRange(data, &path, &opath, cell, &xMin, &xMax, &yMin, &yMax, &xBin, &yBin, &nPlates, &xLow, &yLow);
    
    gStyle->SetOptStat(0);
    TH2::AddDirectory(false);
    TH2F* hm[nPlates];

    // TString outputPath = TString::Format("%s/%i", path.Data(), cell);
    if (!std::filesystem::exists(opath.Data())) std::filesystem::create_directory(opath.Data());
    std::cout << "Combination " << combination << std::endl;
    
    TH2F* hComb = stackHist(data, combination, cell, &hm[0]);
 
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    c->SetGrid();
    hComb->Smooth();
    hComb->Draw("colz");
    c->Update();
    c->Print(Form("%s/sh_%i.gif+180", opath.Data(), combination));
    TObjArray peaks;
    TObjArray txt;
    int ranks[ntag];
    get_peaks(*hComb,peaks,txt,ntag,ranks,bkg);
    hComb->GetZaxis()->SetRangeUser(bkg, hComb->GetMaximum());
    c->Update();
    c->Print(Form("%s/sh_%i.gif+180", opath.Data(), combination));
    drawEllipse(peaks,txt, kBlack);
    c->Update();
    c->Print(Form("%s/sh_%i.gif+180", opath.Data(), combination));

    TH1F *h_long[ntag];
    for(int i=0; i<ntag; i++) {
        h_long[i] = new TH1F(Form("h_long_%i", i+1),Form("Cluster %i;plate;%%segments", i+1), 60, 1, 61);
    }

    for(int p=1; p<=nPlates; p++) { 
        printf("Plate %i\n", p);
        hm[p-1]->Smooth();
        hm[p-1]->Draw("colz");
        drawEllipse(peaks,txt, kBlack);
        count_bins(hm[p-1], peaks, p, &h_long[0], bkg);
        //hm[p-1]->GetZaxis()->SetRangeUser(static_cast<int>(std::ceil((double)bkg/nPlates)), hm[p-1]->GetMaximum());
        c->Update();
        c->Print(Form("%s/sh_%i.gif+12", opath.Data(), combination));
    }
    c->Print(Form("%s/sh_%i.gif++", opath.Data(), combination));

    makeNtuple(combination, cell, &h_long[0], peaks,ranks);

    std::cout << "Time: " << stopWatch.RealTime() << std::endl;
    printMemoryInfo();

    return 0;
}
