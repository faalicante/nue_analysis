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
#include "TMath.h"
#include "TLegend.h"
#include <iostream>
#include <filesystem>
#include <map>
#include <vector>
#include <stdexcept>
#include <cmath>

void printMemoryInfo() {
    ProcInfo_t procInfo;
    gSystem->GetProcInfo(&procInfo);
    Long64_t mem = procInfo.fMemResident;
    std::cout << "Memory used: " << mem/1024 << " MB" << std::endl;
}

// paths
const char* lab = "Napoli";
const int run = 1;
const int brick = 21;

// Parameters
bool print = true;
const int binSize    = 50;   // (um)
const int radius     = 200;  // (um)
const int ntag = 10;
int xMin, xMax, yMin, yMax, xBins, yBins, xLow, yLow;
const int nPlates = 57;
int range;
float bkg = 0;
TString path;
TString opath;
TString ppath;
TString histName;

void getPath(int data, TString* path, TString* opath, TString *ppath, int cell, int* xLow, int* yLow, int* range) {
    if (data == 0) { // Muon simulation
        // *path = "/Users/fabioali/cernbox/shift/muon";
        // *opath = *path;
        *path = TString::Format("/eos/experiment/sndlhc/users/dancc/FEDRA/muon_Euniform_RUN1_FLUKA25/b%06i/cell_reco", brick);
        *opath = "/eos/experiment/sndlhc/users/falicant/RNN/bkg";
        *ppath = TString::Format("%s/%i", opath->Data(), cell);
        *xLow = cell % 18 + 1;
        *yLow = cell / 18 + 1;
        *range = 500;
    }
    else if (data == 1) { // Nue simulation (cell is event)
        // *path = "/Users/fabioali/cernbox/shift/nue_regen";
        // *opath = *path;
        // *ppath = *opath;
        *path = "/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc/nuecc_muon_Euniform_RUN1_FLUKA25/b000021";
        *opath = "/eos/experiment/sndlhc/users/falicant/RNN/signal";
        *ppath = TString::Format("%s/%i", opath->Data(), cell);
        *range = 1000;
    }
    else if (data == 2) { // Real data
        // *path = "/Users/fabioali/cernbox/shift/b121";
        // *opath = *path;
        // *ppath = *opath;
        *path = TString::Format("/eos/experiment/sndlhc/emulsionData/emureco_%s/RUN%i/b%06i/cells", lab, run, brick);
        *opath = TString::Format("/eos/experiment/sndlhc/users/falicant/RUN%i/b%i/shift", run, brick);
        *ppath = TString::Format("/eos/user/f/falicant/RUN%i/brick%i/shifts", run, brick);
        *xLow = cell % 18 + 1;
        *yLow = cell / 18 + 1;
        *range = 4000;
    }
}

TH2F* setRanges1(int data, int cell, TFile** f, int* xMin, int* xMax, int* yMin, int* yMax, int* xBins, int* yBins) {
    TH2F* h2 = (TH2F*)((*f)->Get("XYseg"));
    int fax, fay, lax, lay;
    if (data==1) {
        int nbinsX = h2->GetNbinsX();
        int nbinsY = h2->GetNbinsY();
        fax = nbinsX / 2;
        fay = nbinsY / 2;
        lax = fax;
        lay = fay;
    }
    else {
        fax = h2->FindFirstBinAbove(0,1);
        fay = h2->FindFirstBinAbove(0,2);
        lax = h2->FindLastBinAbove(0,1);
        lay = h2->FindLastBinAbove(0,2);
    }

    *xMin = (int)(h2->GetXaxis()->GetBinLowEdge(fax)) - range;
    *xMax = (int)(h2->GetXaxis()->GetBinUpEdge(lax)) + range;
    *yMin = (int)(h2->GetYaxis()->GetBinLowEdge(fay)) - range;
    *yMax = (int)(h2->GetYaxis()->GetBinUpEdge(lay)) + range;
    *xBins = int((*xMax - *xMin) / binSize);
    *yBins = int((*yMax - *yMin) / binSize);
    h2->GetXaxis()->SetRangeUser(*xMin, *xMax);
    h2->GetYaxis()->SetRangeUser(*yMin, *yMax);
    return h2;
}

void setRanges2(TH2F **hm, float x0, float y0) {
    for (int i = 0; i < nPlates; ++i) {
        hm[i]->GetXaxis()->SetRangeUser(x0-range, x0+range);
        hm[i]->GetYaxis()->SetRangeUser(y0-range, y0+range);
    }
    std::cout << "x range = [" << x0-range << ", " << x0+range << "], y range = [" << y0-range << ", " << y0+range << "]" << std::endl;
}

TH3F* loadH3(TFile *f) {
    TH3F *h3 = nullptr;
    if (f) {
        f->GetObject("XYPseg", h3);
        h3->SetDirectory(0);
    }
    return h3;
}

TH1F* drawSpectrum(TH2F *h2) {
    int nBinsX = h2->GetNbinsX();
    int nBinsY = h2->GetNbinsY();
    TH1F* hSpec = new TH1F("hSpec", "Spectrum;rankbin", 200, 0, 1000);
    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            int content = h2->GetBinContent(i, j);
            if (content > 0) hSpec->Fill(content);
        }
    }
    return hSpec;
}

void poisBkg(TH1F* h, float *bkg) {
    float mpv = h->GetBinLowEdge(h->GetMaximumBin());
    *bkg = mpv+5*std::sqrt(mpv);
    std::cout << "Poisson background: " << *bkg << std::endl;
}

void openFiles(int data, int cell, TFile** f, TH3F** H3cell) {
    TString fileName;
    if (data==1) fileName = TString::Format("%s/b000021.0.0.%i.trk.root", path.Data(), cell+1);
    if (data==0) fileName = TString::Format("%s/cell_%i0_%i0/b%06i/b%06i.0.%i.%i.trk.root", path.Data(), xLow, yLow, brick, brick, xLow, yLow);
    // std::cout << fileName << std::endl;
    *f = TFile::Open(fileName);
    TH2F* H2cell = setRanges1(data, cell, f, &xMin, &xMax, &yMin, &yMax, &xBins, &yBins);
    *H3cell = loadH3(*f);
    H2cell->Smooth();
    TH1F* hSpec2 = drawSpectrum(H2cell);
    poisBkg(hSpec2, &bkg);
}

// void openFiles(int data, int cell, TFile* f[9], TH3F* H3cells[9]) {
//     int idx = 0;
//     for (int yCell = yLow-1; yCell <= yLow+1; yCell ++) {
//         for (int xCell = xLow-1; xCell <= xLow+1; xCell++) {
//             if (xCell < 1 || xCell > 18 || yCell < 1 || yCell > 18) {
//                 f[idx] = nullptr;
//                 H3cells[idx] = nullptr;
//             }
//             else {
//                 TString histFile = TString::Format("%s/cell_%i0_%i0/b%06i/b%06i.0.%i.%i.trk.root", path.Data(), xCell, yCell, brick, brick, xCell, yCell);
//                 // std::cout << histFile << std::endl;
//                 f[idx] = TFile::Open(histFile);
//                 H3cells[idx] = loadH3(f[idx]);
//             }
//             idx++;   
//         }
//     }
//     TH2F* H2cell = setRanges1(data, cell, &f[4], &xMin, &xMax, &yMin, &yMax, &xBins, &yBins);
//     H2cell->Smooth();
//     TH1F* hSpec2 = drawSpectrum(H2cell);
//     poisBkg(hSpec2, &bkg);
// }

TH2F* projectHist(TH3F* h3, int plate) {
    h3->GetEntries();
    h3->GetZaxis()->SetRange(plate+1,plate+1);
    TH2F* h2 = (TH2F*)(h3->Project3D("yx"));
    return h2;
}

TH2F* matrixCells(TH3F* h3, int plate, double shiftX, double shiftY) {
    TH2F* hm = new TH2F(histName, histName, xBins, xMin, xMax, yBins, yMin, yMax);
    TH2F* h2 = projectHist(h3, plate);
    for (int xBin = 1; xBin <= h2->GetNbinsX(); ++xBin) {
        double xCenter = h2->GetXaxis()->GetBinCenter(xBin) + shiftX;
        if (xCenter > xMax || xCenter < xMin) continue;
        for (int yBin = 1; yBin <= h2->GetNbinsY(); ++yBin) {
            double yCenter = h2->GetYaxis()->GetBinCenter(yBin) + shiftY;
            if (yCenter > yMax || yCenter < yMin) continue;
            double content = h2->GetBinContent(xBin, yBin);
            int xBinNew = hm->GetXaxis()->FindBin(xCenter);
            int yBinNew = hm->GetYaxis()->FindBin(yCenter);
            hm->SetBinContent(xBinNew, yBinNew, content);
        }
    }
    delete h2;
    return hm;
}

TH2F* matrixCells(TFile* f[9],  TH3F* H3cells[9], int plate, double shiftX, double shiftY) {
    TH2F* hm = new TH2F(histName, histName, xBins, xMin, xMax, yBins, yMin, yMax);
    TH2F* h2;
    for (int i = 0; i < 9; i++) {
        if (f[i] == nullptr) continue;
        h2 = projectHist(H3cells[i], plate);
        for (int xBin = 1; xBin <= h2->GetNbinsX(); ++xBin) {
            double xCenter = h2->GetXaxis()->GetBinCenter(xBin) + shiftX;
            if (xCenter > xMax || xCenter < xMin) continue;
            for (int yBin = 1; yBin <= h2->GetNbinsY(); ++yBin) {
                double yCenter = h2->GetYaxis()->GetBinCenter(yBin) + shiftY;
                if (yCenter > yMax || yCenter < yMin) continue;
                double content = h2->GetBinContent(xBin, yBin);
                if (content <= 0) continue;
                int xBinNew = hm->GetXaxis()->FindBin(xCenter);
                int yBinNew = hm->GetYaxis()->FindBin(yCenter);
                if (content > hm->GetBinContent(xBinNew, yBinNew)) hm->SetBinContent(xBinNew, yBinNew, content);
            }
        }
        delete h2;
    }
    return hm;
}

TH2F* stackHist(int data, int combination, int cell, TH2F **hm, TString *histName, TH3F *H3cell) {
    TH2F* hComb = new TH2F("XYseg", "XYseg", xBins, xMin, xMax, yBins, yMin, yMax);
    for (int layer = 0; layer < nPlates; ++layer) {
        
        int plate = layer + 1;
        // std::cout << "Shifting plate " << plate << std::endl;
        
        *histName = TString::Format("XYseg_%d", plate);
        hm[layer] = matrixCells(H3cell, plate, 0, 0);
        hComb->Add(hm[layer]);
        hm[layer]->Smooth();
    }
    hComb->Smooth();
    return hComb;
}

TH2F* stackHist(int data, int combination, int cell, TH2F **hm, TString *histName, TFile* ff[9], TH3F* H3cells[9]) {
    TH2F* hComb = new TH2F("XYseg", "XYseg", xBins, xMin, xMax, yBins, yMin, yMax);
    for (int layer = 0; layer < nPlates; ++layer) {
        
        int plate = layer + 1;
        // std::cout << "Shifting plate " << plate << std::endl;
        
        *histName = TString::Format("XYseg_%d", plate);
        hm[layer] = matrixCells(&ff[0], &H3cells[0], plate, 0, 0);
        hComb->Add(hm[layer]);
        hm[layer]->Smooth();
    }
    hComb->Smooth();
    return hComb;
}

int getMax(TH2F &h2, TObjArray &peaks, float bkg) {
    int rankbin = h2.GetMaximum();
    if (rankbin > bkg) {
        Int_t MaxBin = h2.GetMaximumBin();
        Int_t ix,iy,iz;
        h2.GetBinXYZ(MaxBin, ix, iy, iz);
        float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
        float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
        TEllipse  *el = new TEllipse(x,y,radius,radius);
        el->SetFillStyle(0);
        peaks.Add(el);
        int r0 = (int)round((double)radius/binSize);
        for(int iix = ix-r0; iix<=ix+r0; iix++) {
            for(int iiy = iy-r0; iiy<=iy+r0; iiy++) {
                double dx = iix - ix;
                double dy = iiy - iy;
                double distance = (dx*dx + dy*dy)/(r0*r0);
                if (distance <= 1) h2.SetBinContent(iix,iiy,0);
            }
        }
        return rankbin;
    }
    return 0;
}

void get_peaks(TH2F &h2, TObjArray &peaks, int npmax, int *ranks, float bkg) {
    TH2F *h2new = (TH2F*)h2.Clone("get_peaks");
    for(int i=0; i<npmax; i++){
        int rankbin = getMax(*h2new, peaks, bkg);
        ranks[i] = rankbin;
    }
}

double findColScale(TH2F **hm) {
    double maxValue = 0;
    for (int i = 0; i < nPlates; ++i) {
        if (hm[i] == nullptr) continue;
        hm[i]->SetMaximum(-1111);
        double layerMax = hm[i]->GetMaximum();
        double layerMaxBin = hm[i]->GetMaximumBin();
        if (layerMax > maxValue) maxValue = layerMax;
    }
    std::cout << "Max value across all plates: " << maxValue << std::endl;
    return maxValue;
}

void printBW (TH2F **hm, int cell, double zScale, int tag) {
    if (!std::filesystem::exists(ppath.Data())) std::filesystem::create_directory(ppath.Data());
    if (!std::filesystem::exists(TString::Format("%s/%i", ppath.Data(), tag).Data())) {
        std::filesystem::create_directory(TString::Format("%s/%i", ppath.Data(), tag).Data());
    }
    TString imgName;
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    for(int p=1; p<=nPlates; p++) { 
        hm[p-1]->GetZaxis()->SetRangeUser(bkg/nPlates,zScale);
        hm[p-1]->Draw("col0");
        c->Update();
        if (tag>0) imgName = TString::Format("%s/%i/%i_%i_%i.png", ppath.Data(), tag, cell, tag, p);
        else imgName = TString::Format("%s/%i_%i.png", ppath.Data(), cell, p);
        c->Print(imgName);
        c->Clear();
    }
    delete c;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <data>" << argv[1] << " <cell>" << argv[2] << std::endl;
        // data = {0: muon, 1: nue, 2: data}
        // cell for nue is event
    }
    int data = std::atoi(argv[1]);
    int cell = std::atoi(argv[2]);
    
    TStopwatch stopWatch;
    stopWatch.Start();

    gErrorIgnoreLevel = kWarning;
    gROOT->SetBatch(!print);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(52);
    gStyle->SetTitleSize(0, "XYZ");
    gStyle->SetLabelSize(0, "XYZ");
    gStyle->SetFrameLineWidth(0);
    gStyle->SetPadLeftMargin(0);
    gStyle->SetPadRightMargin(0);
    gStyle->SetPadTopMargin(0);
    gStyle->SetPadBottomMargin(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gROOT->ForceStyle();
    
    getPath(data, &path, &opath, &ppath, cell, &xLow, &yLow, &range);
    
    TFile *f, *ff[9];
    TH3F *H3cell, *H3cells[9];

    if (data == 0 || data == 1) {
        openFiles(data, cell, &f, &H3cell);
    }
    // else {
    //     openFiles(data, cell, &ff[0], &H3cells[0]);
    // }

    TH2F* hComb;
    TH2F* hProc;
    TH2F *hm[nPlates];
    TH2::AddDirectory(false);
    int combination=1300;
    stopWatch.Continue();
    
    if (data == 0 || data == 1) hComb = stackHist(data, combination, cell, &hm[0], &histName, H3cell);
    else hComb = stackHist(data, combination, cell, &hm[0], &histName, &ff[0], &H3cells[0]);
    TH1F *hSpec2, *hSpec3;
    TObjArray peaks;
    int ranks[ntag];
    get_peaks(*hComb,peaks,ntag,ranks,bkg);
    
    
    if (data==0) {
        int np = peaks.GetEntries();
        for(int j=0; j<np; j++) {
            TEllipse *el = ((TEllipse*)(peaks.At(j)));
            float x0 = el->GetX1();
            float y0 = el->GetY1();
            std::cout << "Peak " << j+1 << ": x = " << x0 << ", y = " << y0 << std::endl;
            setRanges2(&hm[0], x0, y0);
            double zScale = findColScale(&hm[0]);
            printBW(&hm[0], cell, zScale, j+1);
        }
    }
    
    if (data==1) {
        // double zScale = findColScale(&hm[0]);
        printBW(&hm[0], cell, findColScale(&hm[0]), 0);
    }

    // if (!std::filesystem::exists(ppath.Data())) std::filesystem::create_directory(ppath.Data());
    
    // TCanvas *c = new TCanvas("c", "c", 800, 800);
    // for(int p=1; p<=nPlates; p++) { 
    //     hm[p-1]->Draw("col0");
    //     hm[p-1]->GetZaxis()->SetRangeUser(bkg/nPlates,zScale);
    //     c->Update();
    //     c->Print(Form("%s/%i_%i.png", ppath.Data(), cell, p));
    //     c->Clear();
    // }

    delete hComb;
    
    for(int p=1; p<=nPlates; p++) { 
        delete hm[p-1];
    }
    std::cout << "---------------------" << std::endl;
    
    if (data == 1 || data == 0) {
        f->Close();
        delete H3cell;
    }
    else {
        for (int i = 0; i < 9; i++) {
            if (ff[i] != nullptr) ff[i]->Close();
            delete H3cells[i];
        }
    }

    std::cout << "Time: " << round(stopWatch.RealTime()) << std::endl;
    printMemoryInfo();
    
    return 0;
}
