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
    ProcInfo_t procInfo;
    gSystem->GetProcInfo(&procInfo);
    Long64_t mem = procInfo.fMemResident;
    std::cout << "Memory used: " << mem/1024 << " MB" << std::endl;
}


// Parameters
bool print = false;
const int binSize    = 50;   // (um)
const int shiftRange = 50;   // (mrad)
const int shiftStep  = 2;    // (mrad)
const int radius     = 300;  // (um)
const int ntag = 20;
int xMin, xMax, yMin, yMax, xBins, yBins, xLow, yLow;
int nPlates, stepZ;
float bkg = 0;
TString path;
TString opath;
TString histName;

void getPath(int data, TString* path, TString* opath, int cell, int* nPlates, int* xLow, int* yLow, int* stepZ) {
    if (data == 0) { // Muon simulation
        // *path = "/Users/fabioali/cernbox/shift/muon";
        // *opath = *path;
        *path = "/eos/experiment/sndlhc/MonteCarlo/FEDRA/muon1.3E5/cell_reco";
        *opath = TString::Format("/eos/experiment/sndlhc/users/falicant/shift_muon/%i", cell);
        *nPlates = 60;
        *stepZ = 1315;
        *xLow = cell % 21 + 1;
        *yLow = cell / 21 + 1;
    }
    else if (data == 1) { // Nue simulation (cell is event)
        // *path = "/Users/fabioali/cernbox/shift/nue_muon";
        // *opath = *path;
        *path = "/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc/nuecc_muon1.3E5/b000021";
        *opath = TString::Format("/eos/experiment/sndlhc/users/falicant/shift_nue/%i", cell);
        *nPlates = 60;
        *stepZ = 1315;
    }
    else if (data == 2) { // Real data
        // *path = "/Users/fabioali/cernbox/shift/b121";
        // *opath = *path;
        *path = "/eos/experiment/sndlhc/emulsionData/emureco_Napoli/RUN1/b000121/cells";
        *opath = TString::Format("/eos/experiment/sndlhc/users/falicant/RUN1/b121/shift/%i", cell);
        *nPlates = 57;
        *stepZ = 1350;
        *xLow = cell % 18 + 1;
        *yLow = cell / 18 + 1;
    }
}

void setRanges(int cell, TFile** f, int* xMin, int* xMax, int* yMin, int* yMax, int* xBins, int* yBins, float* bkg) {
    TH2F* h2 = (TH2F*)((*f)->Get("XYseg"));
    int fax = h2->FindFirstBinAbove(0,1);
    int fay = h2->FindFirstBinAbove(0,2);
    int lax = h2->FindLastBinAbove(0,1);
    int lay = h2->FindLastBinAbove(0,2);
    *xMin = (int)(h2->GetXaxis()->GetBinLowEdge(fax)) - 1000;
    *xMax = (int)(h2->GetXaxis()->GetBinUpEdge(lax)) + 1000;
    *yMin = (int)(h2->GetYaxis()->GetBinLowEdge(fay)) - 1000;
    *yMax = (int)(h2->GetYaxis()->GetBinUpEdge(lay)) + 1000;
    *xBins = int((*xMax - *xMin) / binSize);
    *yBins = int((*yMax - *yMin) / binSize);
    *bkg = 1.5 * h2->Integral(fax, lax, fay, lay) / (lax - fax + 1) / (lay - fay + 1);
    std::cout << "Average background " << *bkg << std::endl;
    std::cout << *xMin << " " << *xMax << " " << *yMin << " " << *yMax << std::endl;
    std::cout << *xBins << " " << *yBins << std::endl;
    delete h2;
}

void openFiles(int cell, TFile** f) {
    TString fileName = TString::Format("%s/b000021.0.0.%i.trk.root", path.Data(), cell+1);
    std::cout << fileName << std::endl;
    *f = TFile::Open(fileName);
    setRanges(cell, f, &xMin, &xMax, &yMin, &yMax, &xBins, &yBins, &bkg);
}

void openFiles(int data, int cell, TFile* f[9]) {
    int idx = 0;
    if (data == 0) {
        for (int yCell = yLow-1; yCell <= yLow+1; yCell ++) {
            for (int xCell = xLow-1; xCell <= xLow+1; xCell++) {
                if (xCell < 1 || xCell > 19 || yCell < 1 || yCell > 19) f[idx] = nullptr;
                else {
                    int iCell = (yCell-1)*21 + (xCell-1);
                    TString histFile = TString::Format("%s/%i/b000021/b000021.0.0.0.trk.root", path.Data(), iCell);
                    std::cout << histFile << std::endl;
                    f[idx] = TFile::Open(histFile);
                }
                idx++;   
            }
        }
    }
    else if (data == 2) {
        for (int yCell = yLow-1; yCell <= yLow+1; yCell ++) {
            for (int xCell = xLow-1; xCell <= xLow+1; xCell++) {
                if (xCell < 1 || xCell > 18 || yCell < 1 || yCell > 18) f[idx] = nullptr;
                else {
                    TString histFile = TString::Format("%s/cell_%i0_%i0/b000021/b000021.0.0.0.trk.root", path.Data(), xCell, yCell);
                    std::cout << histFile << std::endl;
                    f[idx] = TFile::Open(histFile);
                }
                idx++;   
            }
        }
    }
    setRanges(cell, &f[4], &xMin, &xMax, &yMin, &yMax, &xBins, &yBins, &bkg);
}

TH2F* projectHist(TFile* f, int plate) {
    TH3F* h3 = (TH3F*)(f->Get("XYPseg"));
    h3->GetZaxis()->SetRange(plate+1,plate+1);
    TH2F* h2 = (TH2F*)(h3->Project3D("yx"));
    delete h3;
    return h2;
}

TH2F* matrixCells(TFile* f, int plate, double shiftX, double shiftY) {
    TH2F* hm = new TH2F(histName, histName, xBins, xMin, xMax, yBins, yMin, yMax);
    TH2F* h2 = projectHist(f, plate);
    for (int xBin = 1; xBin <= h2->GetNbinsX(); ++xBin) {
        double xCenter = h2->GetXaxis()->GetBinCenter(xBin) + shiftX;
        if (xCenter <= xMax && xCenter >= xMin) {
            for (int yBin = 1; yBin <= h2->GetNbinsY(); ++yBin) {
                double yCenter = h2->GetYaxis()->GetBinCenter(yBin) + shiftY;
                if (yCenter <= yMax && yCenter >= yMin) {
                    double content = h2->GetBinContent(xBin, yBin);
                    int xBinNew = hm->GetXaxis()->FindBin(xCenter);
                    int yBinNew = hm->GetYaxis()->FindBin(yCenter);
                    hm->SetBinContent(xBinNew, yBinNew, content);
                }
            }
        }
    }
    return hm;
}

TH2F* matrixCells(TFile* f[9], int plate, double shiftX, double shiftY) {
    TH2F* hm = new TH2F(histName, histName, xBins, xMin, xMax, yBins, yMin, yMax);
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

TH2F* stackHist(int data, int combination, int cell, TH2F **hm, TString *histName) {
    TFile *f, *ff[9];

    if (data == 1) {
        openFiles(cell, &f);
    }
    else {
        openFiles(data, cell, &ff[0]);
    }

    double shiftTX = (combination % (shiftRange+1)) * shiftStep - shiftRange;
    double shiftTY = (combination / (shiftRange+1)) * shiftStep - shiftRange;
    // combination = (shiftTY + shiftRange)/shiftStep * (shiftRange + 1) + (shiftTX + shiftRange)/shiftStep
    std::cout << "Shift TX: " << shiftTX << " mrad, Shift TY: " << shiftTY << " mrad" << std::endl;
    
    TH2F* hComb = new TH2F("XYseg", "XYseg", xBins, xMin, xMax, yBins, yMin, yMax);
    for (int layer = 0; layer < nPlates; ++layer) {
        int plate = layer + 1;
        // std::cout << "Shifting plate " << plate << std::endl;

        double shiftX = shiftTX / 1000.0 * stepZ * layer;
        double shiftY = shiftTY / 1000.0 * stepZ * layer;
        
        *histName = TString::Format("XYseg_%d", plate);
        if (data == 1) hm[layer] = matrixCells(f, plate, shiftX, shiftY);
        else hm[layer] = matrixCells(&ff[0], plate, shiftX, shiftY);
        hComb->Add(hm[layer]);
    }
    if (data == 1) f->Close();
    else {
        for (int i = 0; i < 9; i++) {
            if (ff[i] != nullptr) ff[i]->Close();
        }
    }
    return hComb;
}

int getMax(TH2F &h2, TObjArray &peaks, TObjArray &txt, float bkg) {
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
        TText  *t = new TText(x,y+300,Form("%d",peaks.GetEntries()));
        t->SetTextSize(0.03);
        txt.Add(t);
        int r0 = (int)round((double)radius/binSize);
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

void get_peaks(TH2F &h2, TObjArray &peaks, TObjArray &txt, int npmax, int *ranks, float bkg) {
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
        ((TText*)(txt.At(j)))->Draw();
    }
}

void getEntriesInEllipse(TH2F &h2, TEllipse &el, int *entries, float bkg) {
    int nBinsX = h2.GetNbinsX();
    int nBinsY = h2.GetNbinsY();
    float x0 = el.GetX1();
    float y0 = el.GetY1();
    float r1 = el.GetR1();
    *entries = 0;
    float bkg_p = static_cast<int>(std::ceil((double)bkg/nPlates));
    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            float x = h2.GetXaxis()->GetBinCenter(i);
            float y = h2.GetYaxis()->GetBinCenter(j);
            float distance = (pow(x-x0,2)+pow(y-y0,2))/pow(r1,2);
            if (distance <= 1) {
                int binContent = h2.GetBinContent(i, j);
                if(binContent > bkg_p) {
                    *entries = *entries + binContent - bkg_p;
                }
            }
        } 
    }
}

void count_bins(TH2F *h2, TObjArray &peaks, int plate, TH1F **h_long, float bkg) {
    int np = peaks.GetEntries();
    int entries;

    for(int i=0; i<np; i++) {
        TEllipse *el = ((TEllipse*)(peaks.At(i)));
        getEntriesInEllipse(*h2, *el, &entries, bkg);
        h_long[i]->SetBinContent(plate, entries);
    }
}

void makePlots(int cell, int combination, TCanvas *c, int np, int npmax, TH1F *h_long, int *maxPeak, int *maxPlate) {
    int idx = (np%3) +1;
    if (idx==1) c->Clear("D");
    c->cd(idx)->SetGrid(1,0);
    *maxPeak = h_long->GetMaximum();
    *maxPlate = h_long->GetMaximumBin();
    h_long->SetLineColor(1);
    h_long->SetLineWidth(2);
    h_long->Draw("hist");
    if (print) {
        if(np == 2)                 c->Print(Form("%s/longitudinal_xz_%i_%i.pdf(", opath.Data(), cell, combination), "pdf");
        else if(np == npmax-1)      c->Print(Form("%s/longitudinal_xz_%i_%i.pdf)", opath.Data(), cell, combination), "pdf");
        else if(idx == 3 && np >2 ) c->Print(Form("%s/longitudinal_xz_%i_%i.pdf" , opath.Data(), cell, combination), "pdf");
    }
}

void findStart(TH1F* h_long, int *firstPlate, int *lastPlate) { // problem with holes
    float entries = h_long->Integral();
    *firstPlate = h_long->FindFirstBinAbove(entries * 0.005);
    *lastPlate = h_long->FindLastBinAbove(entries * 0.005);
    // std::cout << entries << " " << *firstPlate << std::endl;
}

void makeNtuple(int combination, int cell, TH1F **h_long, TObjArray &peaks, int *ranks) {
    TString outputFileName = TString::Format("%s/peaks_%i_%i.root", opath.Data(), cell, combination);
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
        makePlots(cell, combination, c2, i, np, h_long[i], &maxPeak, &maxPlate);
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
        // cell for nue is event
        return 1;
    }
    int data = std::atoi(argv[1]);
    int combination = std::atoi(argv[2]);
    int cell = std::atoi(argv[3]);

    TStopwatch stopWatch;
    stopWatch.Start();

    getPath(data, &path, &opath, cell, &nPlates, &xLow, &yLow, &stepZ);
    
    gStyle->SetOptStat(0);
    TH2::AddDirectory(false);
    TH2F* hm[nPlates];

    if (!std::filesystem::exists(opath.Data())) std::filesystem::create_directory(opath.Data());
    std::cout << "Combination " << combination << std::endl;
    
    TH2F* hComb = stackHist(data, combination, cell, &hm[0], &histName);
 
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    c->SetGrid();
    hComb->Smooth();
    hComb->Draw("colz");
    c->Update();
    if (print) c->Print(Form("%s/sh_%i_%i.gif+180", opath.Data(), cell, combination));
    TObjArray peaks;
    TObjArray txt;
    int ranks[ntag];
    get_peaks(*hComb,peaks,txt,ntag,ranks,bkg);
    hComb->GetZaxis()->SetRangeUser(bkg, hComb->GetMaximum());
    c->Update();
    if (print) c->Print(Form("%s/sh_%i_%i.gif+180", opath.Data(), cell, combination));
    drawEllipse(peaks,txt, kBlack);
    c->Update();
    if (print) c->Print(Form("%s/sh_%i_%i.gif+180", opath.Data(), cell, combination));

    TH1F *h_long[ntag];
    for(int i=0; i<ntag; i++) {
        h_long[i] = new TH1F(Form("h_long_%i", i+1),Form("Cluster %i;plate;%%segments", i+1), 60, 1, 61);
    }
    
    for(int p=1; p<=nPlates; p++) { 
        // printf("Tagging plate %i\n", p);
        hm[p-1]->Smooth();
        hm[p-1]->Draw("colz");
        drawEllipse(peaks,txt, kBlack);
        count_bins(hm[p-1], peaks, p, &h_long[0], bkg);
        // hm[p-1]->GetZaxis()->SetRangeUser(static_cast<int>(std::ceil((double)bkg/nPlates)), hm[p-1]->GetMaximum());
        c->Update();
        if (print) c->Print(Form("%s/sh_%i_%i.gif+12", opath.Data(), cell, combination));
    }
    if (print) c->Print(Form("%s/sh_%i_%i.gif++", opath.Data(), cell, combination));

    makeNtuple(combination, cell, &h_long[0], peaks,ranks);

    std::cout << "Time: " << stopWatch.RealTime() << std::endl;
    printMemoryInfo();

    return 0;
}