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
#include "TSpectrum.h"
#include "TMath.h"
#include "TLegend.h"
#include "TF1.h"
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
const int brick = 121;

// Parameters
bool print = false;
const int binSize    = 50;   // (um)
const int shiftRange = 50;   // (mrad) //enlarge with bigger step
const int shiftStep  = 2;    // (mrad)
const int radius     = 200;  // (um)
const int ntag = 100;
int xMin, xMax, yMin, yMax, xBins, yBins, xLow, yLow;
const int nPlates = 57;
const int stepZ = 1350;
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
        *path = "/eos/experiment/sndlhc/users/dancc/FEDRA/muon_regenRUN1/cell_reco";
        *opath = TString::Format("/eos/experiment/sndlhc/users/falicant/shift_muon/%i", cell);
        *ppath = "/eos/user/f/falicant/shift/muon_regen";
        *xLow = cell % 18 + 1;
        *yLow = cell / 18 + 1;
        *range = 4000;
    }
    else if (data == 1) { // Nue simulation (cell is event)
        // *path = "/Users/fabioali/cernbox/shift/nue_regen";
        // *opath = *path;
        // *ppath = *opath;
        *path = "/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc/nuecc_muon_regenRUN1/b000021";
        *opath = "/eos/experiment/sndlhc/users/falicant/shift_nue_regen_100";
        *ppath = "/eos/user/f/falicant/shift/nue_regen";
        *range = 0;
    }
    else if (data == 2) { // Real data
        // *path = "/Users/fabioali/cernbox/shift/b121";
        // *opath = *path;
        // *ppath = *opath;
        *path = TString::Format("/eos/experiment/sndlhc/emulsionData/emureco_%s/RUN%i/b%06i/cells", lab, run, brick);
        *opath = TString::Format("/eos/experiment/sndlhc/users/falicant/RUN%i/b%i/shift", run, brick, cell);
        *ppath = TString::Format("/eos/user/f/falicant/RUN%i/brick%i/shifts", run, brick);
        *xLow = cell % 18 + 1;
        *yLow = cell / 18 + 1;
        *range = 4000;
    }
}

void setRanges(int cell, TFile** f, int* xMin, int* xMax, int* yMin, int* yMax, int* xBins, int* yBins) {
    TH2F* h2 = (TH2F*)((*f)->Get("XYseg"));
    int fax = h2->FindFirstBinAbove(0,1);
    int fay = h2->FindFirstBinAbove(0,2);
    int lax = h2->FindLastBinAbove(0,1);
    int lay = h2->FindLastBinAbove(0,2);
    *xMin = (int)(h2->GetXaxis()->GetBinLowEdge(fax)) - range;
    *xMax = (int)(h2->GetXaxis()->GetBinUpEdge(lax)) + range;
    *yMin = (int)(h2->GetYaxis()->GetBinLowEdge(fay)) - range;
    *yMax = (int)(h2->GetYaxis()->GetBinUpEdge(lay)) + range;
    *xBins = int((*xMax - *xMin) / binSize);
    *yBins = int((*yMax - *yMin) / binSize);
    delete h2;
}

TH3F* loadH3(TFile *f) {
    TH3F *h3 = nullptr;
    if (f) {
        f->GetObject("XYPseg", h3);
        h3->SetDirectory(0);
    }
    return h3;
}

TH2F* loadH2(TFile *f) {
    TH2F *h2 = nullptr;
    if (f) {
        f->GetObject("XYseg", h2);
        h2->SetDirectory(0);
    }
    return h2;
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

TF1* fitBackground(TH1F* h, float *bkg) {
    float maxBin = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
    TF1 *f = new TF1("f", "gaus(0)+ [3]*exp(-[4]*(x-[1]))/(1+exp(-[5]*(x-[1])))", maxBin-200,maxBin+500); 
    f->SetParameters(
        h->Integral(),     // Gaussian amplitude
        maxBin,            // peak position
        50,                // sigma
        h->Integral(),     // exponential amplitudes
        0.005              // decay constant
    );
    f->SetParLimits(2, 1, 500);   // sigma > 0
    f->SetParLimits(4, 1e-5, 1);  // lambda > 0
    h->Fit(f, "RMQ");
    *bkg = f->GetParameter(1)+5*f->GetParameter(2);
    std::cout << "Fitted background: " << *bkg << std::endl;
    return f;
}

void openFiles(int cell, TFile** f, TH3F** H3cell) {
    TString fileName = TString::Format("%s/b000021.0.0.%i.trk.root", path.Data(), cell+1);
    // std::cout << fileName << std::endl;
    *f = TFile::Open(fileName);
    *H3cell = loadH3(*f);
    TH2F* H2cell = loadH2(*f);
    H2cell->Smooth();
    TH1F* hSpec2 = drawSpectrum(H2cell);
    TF1* fit = fitBackground(hSpec2, &bkg);
    setRanges(cell, f, &xMin, &xMax, &yMin, &yMax, &xBins, &yBins);
}

void openFiles(int data, int cell, TFile* f[9], TH3F* H3cells[9]) {
    int idx = 0;
    for (int yCell = yLow-1; yCell <= yLow+1; yCell ++) {
        for (int xCell = xLow-1; xCell <= xLow+1; xCell++) {
            if (xCell < 1 || xCell > 18 || yCell < 1 || yCell > 18) {
                f[idx] = nullptr;
                H3cells[idx] = nullptr;
            }
            else {
                TString histFile = TString::Format("%s/cell_%i0_%i0/b000021/b000021.0.0.0.trk.root", path.Data(), xCell, yCell);
                // std::cout << histFile << std::endl;
                f[idx] = TFile::Open(histFile);
                H3cells[idx] = loadH3(f[idx]);
            }
            idx++;   
        }
    }
    TH2F* H2cell = loadH2(f[4]);
    H2cell->Smooth();
    TH1F* hSpec2 = drawSpectrum(H2cell);
    TF1* fit = fitBackground(hSpec2, &bkg);
    setRanges(cell, &f[4], &xMin, &xMax, &yMin, &yMax, &xBins, &yBins);
}

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
    double shiftTX = (combination % (shiftRange+1)) * shiftStep - shiftRange;
    double shiftTY = (combination / (shiftRange+1)) * shiftStep - shiftRange;
    // combination = (shiftTY + shiftRange)/shiftStep * (shiftRange + 1) + (shiftTX + shiftRange)/shiftStep
    // std::cout << "Shift TX: " << shiftTX << " mrad, Shift TY: " << shiftTY << " mrad" << std::endl;
    
    TH2F* hComb = new TH2F("XYseg", "XYseg", xBins, xMin, xMax, yBins, yMin, yMax);
    for (int layer = 0; layer < nPlates; ++layer) {
        
        int plate = layer + 1;
        // std::cout << "Shifting plate " << plate << std::endl;
        
        double shiftX = shiftTX / 1000.0 * stepZ * layer;
        double shiftY = shiftTY / 1000.0 * stepZ * layer;
        
        *histName = TString::Format("XYseg_%d", plate);
        hm[layer] = matrixCells(H3cell, plate, shiftX, shiftY);
        hComb->Add(hm[layer]);
    }
    return hComb;
}

TH2F* stackHist(int data, int combination, int cell, TH2F **hm, TString *histName, TFile* ff[9], TH3F* H3cells[9]) {
    double shiftTX = (combination % (shiftRange+1)) * shiftStep - shiftRange;
    double shiftTY = (combination / (shiftRange+1)) * shiftStep - shiftRange;
    // combination = (shiftTY + shiftRange)/shiftStep * (shiftRange + 1) + (shiftTX + shiftRange)/shiftStep
    // std::cout << "Shift TX: " << shiftTX << " mrad, Shift TY: " << shiftTY << " mrad" << std::endl;
    
    TH2F* hComb = new TH2F("XYseg", "XYseg", xBins, xMin, xMax, yBins, yMin, yMax);
    for (int layer = 0; layer < nPlates; ++layer) {
        
        int plate = layer + 1;
        // std::cout << "Shifting plate " << plate << std::endl;
        
        double shiftX = shiftTX / 1000.0 * stepZ * layer;
        double shiftY = shiftTY / 1000.0 * stepZ * layer;
        
        *histName = TString::Format("XYseg_%d", plate);
        hm[layer] = matrixCells(&ff[0], &H3cells[0], plate, shiftX, shiftY);
        hComb->Add(hm[layer]);
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
        t->SetTextSize(0.02);
        txt.Add(t);
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

TH2F* get_peaks(TH2F &h2, TObjArray &peaks, TObjArray &txt, int npmax, int *ranks, float bkg) {
    TH2F *h2new = (TH2F*)h2.Clone("get_peaks");
    for(int i=0; i<npmax; i++){
        int rankbin = getMax(*h2new, peaks, txt, bkg);
        ranks[i] = rankbin;
    }
    return h2new;
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
            double dx = x - x0;
            double dy = y - y0;
            double distance = (dx*dx + dy*dy)/(r1*r1);
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

void makePlots(int cell, int combination, TCanvas *c, int np, int npmax, TH1F *h_long) {
    int idx = (np%3) +1;
    if (idx==1) c->Clear("D");
    c->cd(idx)->SetGrid(1,0);
    h_long->SetLineColor(1);
    h_long->SetLineWidth(2);
    h_long->Draw("hist");
    if (print) {
        if(np == 2)                 c->Print(Form("%s/longitudinal_xz_%i_%i.pdf(", ppath.Data(), cell, combination), "pdf");
        else if(np == npmax-1)      c->Print(Form("%s/longitudinal_xz_%i_%i.pdf)", ppath.Data(), cell, combination), "pdf");
        else if(idx == 3 && np >2 ) c->Print(Form("%s/longitudinal_xz_%i_%i.pdf" , ppath.Data(), cell, combination), "pdf");
    }
}

void findStart(TH1F* h_long, int *firstPlate, int *lastPlate, int *nfound, int *maxPeak, int *maxPlate) {
    float entries = h_long->GetMaximum();
    *maxPeak = h_long->GetMaximum();
    *maxPlate = h_long->GetMaximumBin();
    *firstPlate = h_long->FindFirstBinAbove(entries * 0.1);
    *lastPlate = h_long->FindLastBinAbove(entries * 0.1);
    TSpectrum *s = new TSpectrum(4);
    *nfound = s->Search(h_long, 5, "nobackground", 0.05);
}

void makeNtuple(TFile* outputFile, TNtuple* ntuple, int combination, int cell, TH1F **h_long, TObjArray &peaks, int *ranks) {
    TString canvasName = TString::Format("c2_%i.root", combination);
    TCanvas *c2 = new TCanvas(canvasName, canvasName, 1500, 1500);
    c2->Divide(1,3);

    int np = peaks.GetEntries();
    int entries, firstPlate, lastPlate, nfound, maxPeak, maxPlate;
    
    for(int i=0; i<np; i++) {
        TEllipse *el = ((TEllipse*)(peaks.At(i)));
        float x = el->GetX1();
        float y = el->GetY1();
        if (print) makePlots(cell, combination, c2, i, np, h_long[i]);
        findStart(h_long[i], &firstPlate, &lastPlate, &nfound, &maxPeak, &maxPlate);
        int nseg = h_long[i]->Integral(firstPlate, lastPlate);
        if (print) h_long[i]->Write();
        ntuple->Fill(cell, combination, i+1, x, y, firstPlate, lastPlate, maxPeak, maxPlate, nseg, nfound, ranks[i], bkg);
    }
    delete c2;
}



int main(int argc, char* argv[]) {
    gROOT->SetBatch(!print);
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <data>" << argv[1] << " <cell>" << argv[2] << std::endl;
        // data = {0: muon, 1: nue, 2: data}
        // cell for nue is event
        return 1;
    }
    int data = std::atoi(argv[1]);
    int cell = std::atoi(argv[2]);
    
    TStopwatch stopWatch;
    stopWatch.Start();
    
    getPath(data, &path, &opath, &ppath, cell, &xLow, &yLow, &range);
    
    TFile *f, *ff[9];
    TH3F *H3cell, *H3cells[9];

    if (data == 1) {
        openFiles(cell, &f, &H3cell);
    }
    else {
        openFiles(data, cell, &ff[0], &H3cells[0]);
    }

    gStyle->SetOptStat(0);
    
    TString outputFileName = TString::Format("%s/peaks_%i.root", opath.Data(), cell);
    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    TNtuple *ntuple = new TNtuple("showers","tagged showers","cell:combination:tag:x:y:start:end:peak:maxplate:nseg:nfound:rankbin:bkg");
    
    // if (!std::filesystem::exists(opath.Data())) std::filesystem::create_directory(opath.Data());
    TH2F* hComb;
    TH2F* hProc;
    TH2F *hm[nPlates];
    TH2::AddDirectory(false);
    for(int combination = 0; combination < ((shiftRange+1)*(shiftRange+1)); combination++) {
        if (combination!=1300) continue;
        stopWatch.Continue();
        
        std::cout << "Combination " << combination << std::endl;
        
        if (data == 1) hComb = stackHist(data, combination, cell, &hm[0], &histName, H3cell);
        else hComb = stackHist(data, combination, cell, &hm[0], &histName, &ff[0], &H3cells[0]);
        TString canvasName = TString::Format("c_%i.root", combination);
        TCanvas *c = new TCanvas(canvasName, canvasName, 800, 800);
        TH1F* hSpec1, *hSpec2, *hSpec3;
        if (print) {
            hSpec1 = drawSpectrum(hComb);
            hSpec1->SetLineColor(kBlue);
        }
        hComb->Smooth();
        if (print) {
            hSpec2 = drawSpectrum(hComb);
            hSpec2->SetLineColor(kGreen);
            // TF1* fit = fitBackground(hSpec2, &bkg);
            c->SetGrid();
            hComb->Draw("colz");
            c->Update();
            c->Print(Form("%s/sh_%i_%i.gif+180", ppath.Data(), cell, combination));
        }
        TObjArray peaks;
        TObjArray txt;
        int ranks[ntag];
        hProc = get_peaks(*hComb,peaks,txt,ntag,ranks,bkg);
        if (print) {
            hSpec3 = drawSpectrum(hProc);
            hSpec3->SetLineColor(kRed);
            drawEllipse(peaks,txt, kBlack);
            c->Update();
            c->Print(Form("%s/sh_%i_%i.gif+180", ppath.Data(), cell, combination));
            hComb->GetZaxis()->SetRangeUser(bkg, hComb->GetMaximum());
            c->Update();
            c->Print(Form("%s/sh_%i_%i.gif+180", ppath.Data(), cell, combination));
        }

        TH1F *h_long[ntag];
        for(int i=0; i<ntag; i++) {
            h_long[i] = new TH1F(Form("h_long_%i", i+1),Form("Cluster %i;plate;%%segments", i+1), 60, 1, 61);
        }
        
        for(int p=1; p<=nPlates; p++) { 
            // printf("Tagging plate %i\n", p);
            hm[p-1]->Smooth();
            hm[p-1]->Draw("colz");
            count_bins(hm[p-1], peaks, p, &h_long[0], bkg);
            if (print) {
                drawEllipse(peaks,txt, kBlack);
                // hm[p-1]->GetZaxis()->SetRangeUser(static_cast<int>(std::ceil((double)bkg/nPlates)), hm[p-1]->GetMaximum());
                c->Update();
                c->Print(Form("%s/sh_%i_%i.gif+12", ppath.Data(), cell, combination));
            }
        }
        if (print) c->Print(Form("%s/sh_%i_%i.gif++", ppath.Data(), cell, combination));

        makeNtuple(outputFile, ntuple, combination, cell, &h_long[0], peaks,ranks);

        std::cout << "Time: " << round(stopWatch.RealTime()) << std::endl;
        printMemoryInfo();
        delete hComb;
        delete c;
        for(int p=1; p<=nPlates; p++) { 
            delete hm[p-1];
        }

        if (print) {
            TString canvasNameSp = TString::Format("csp_%i.root", combination);
            TCanvas *cSp = new TCanvas(canvasNameSp, canvasNameSp, 800, 600);
            cSp->cd()->SetLogy();
            cSp->SetGrid();
            hSpec3->Draw("hist ");
            hSpec2->Draw("hist same");
            hSpec1->Draw("hist same");
            // fit->Draw("same"); //TLine for bkg
            TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
            leg->AddEntry(hSpec1, "Original", "l");
            leg->AddEntry(hSpec2, "Smoothed", "l");
            leg->AddEntry(hSpec3, "After peak search", "l");
            leg->Draw();
            cSp->Print(Form("%s/spectrum_%i_%i.pdf", ppath.Data(), cell, combination), "pdf");
            delete cSp;
            delete hSpec1;
            delete hSpec2;
            delete hSpec3;
        }
        std::cout << "---------------------" << std::endl;
    }
    if (data == 1) {
        f->Close();
        delete H3cell;
    }
    else {
        for (int i = 0; i < 9; i++) {
            if (ff[i] != nullptr) ff[i]->Close();
            delete H3cells[i];
        }
    }
    outputFile->Write();
    outputFile->Close();
    return 0;
}
