#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TKey.h"
#include "TROOT.h"
#include "TList.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>
#include <filesystem>
#include <map>
#include <vector>
#include <stdexcept>
#include <cmath>

// Function to load histograms
std::map<std::string, TH2F*> loadHists(const char* histFile) {
    TFile* f = TFile::Open(histFile);
    std::map<std::string, TH2F*> histList;
    TList* keyList = f->GetListOfKeys();
    TIter next(keyList);
    TKey* key;
    TString exclude = "segments";
    while ((key = (TKey*)next())) {
        TString keyName = key->GetName();
        if (keyName == exclude) continue;
        TH2F* hist = (TH2F*)f->Get(key->GetName());
        hist->SetDirectory(gROOT);
        hist->SetName(key->GetName());
        histList[key->GetName()] = hist;
    }
    if (histList.empty()) throw std::runtime_error("ERROR: histList is empty!");
    f->Close();
    return histList;
}

void printMemoryInfo() {
    ProcInfo_t procInfo;  // Declare a ProcInfo_t structure
    gSystem->GetProcInfo(&procInfo);  // Pass its pointer to GetProcInfo

    Long64_t mem = procInfo.fMemResident;  // Access the memory resident field
    std::cout << "Memory used: " << mem/1024 << " MB" << std::endl;
}

// Parameters
const int binSize    = 50;  //um
const int shiftRange = 50; //mrad
const int shiftStep  = 2;   //mrad
const int stepZ      = 1000 + 350; //um
int xMin, xMax, yMin, yMax, xBin, yBin;
int nPlates;
TString path;
// TString path = "/Users/fabioali/cernbox"; //for local test

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

void setRange(int data, TString* path, int cell, int* xMin, int* xMax, int* yMin, int* yMax, int* xBin, int* yBin, int* nPlates) {
    if (data==0) {
        // *path = "/Users/fabioali/cernbox";
        *path = TString::Format("/eos/experiment/sndlhc/users/falicant/Simulations_sndlhc/muon1E5_simsndlhc/b%06i", cell);
        *nPlates = 60;
        *xMin = 284000;
        *xMax = 304000;
        *yMin = 79000;
        *yMax = 99000;
    }
    else if (data==1) {
        *path = TString::Format("/eos/experiment/sndlhc/users/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b%06i", cell);
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
        // *path = "/Users/fabioali/cernbox";
        *path = "/eos/experiment/sndlhc/users/falicant/RUN1/b121/hist";
        *nPlates = 57;
        const int xLow = cell / 18;
        const int yLow = cell % 18;
        *xMin = xLow*10000 + 1000;
        *xMax = xLow*10000 + 19000;
        *yMin = yLow*10000 + 1000;
        *yMax = yLow*10000 + 19000;
    }
    *xBin = int((*xMax - *xMin) / binSize);
    *yBin = int((*yMax - *yMin) / binSize);
}

void openFiles(int cell, TFile* f[9]) {
    int idx = 0;
    for (int yCell = cell-18; yCell <= cell+18; yCell +=18) {
        for (int xCell = yCell-1; xCell <= yCell+1; xCell++) {
            if (xCell < 0) continue;
            TString histFile = TString::Format("%s/hist_XYP_b121_%i.root", path.Data(), xCell);
            f[idx] = TFile::Open(histFile);
            idx++;   
        }
    }
}

TH2F* matrixCells(TFile* f[9], int plate) {
    TString histName = TString::Format("XYseg_%d", plate);
    TString histTitle = TString::Format("XYPseg_%d", plate);
    TList *list = new TList;
    TH2F* hm = new TH2F(histName, histTitle, xBin, xMin, xMax, yBin, yMin, yMax);
    TH2F* h2;
    for (int i = 0; i < 9; i++) {
        h2 = (TH2F*)f[i]->Get(histName);
        list->Add(h2);
        h2->SetDirectory(0);
    }
    hm->Merge(list);
    list->Delete();
    delete list;
    return hm;
}

TH2F* cropHist(TH2F* h2, double shiftX, double shiftY) {
    TH2F* hCrop = new TH2F(h2->GetTitle(), h2->GetTitle(), xBin, xMin, xMax, yBin, yMin, yMax);
    for (int xBin = 1; xBin <= h2->GetNbinsX(); ++xBin) {
        double xCenter = h2->GetXaxis()->GetBinCenter(xBin) + shiftX;
        if (xCenter <= xMax && xCenter >= xMin) {
            for (int yBin = 1; yBin <= h2->GetNbinsY(); ++yBin) {
                double yCenter = h2->GetYaxis()->GetBinCenter(yBin) + shiftY;
                if (yCenter <= yMax && yCenter >= yMin) {
                    double content = h2->GetBinContent(xBin, yBin);
                    int xBinNew = hCrop->GetXaxis()->FindBin(xCenter);
                    int yBinNew = hCrop->GetYaxis()->FindBin(yCenter);
                    hCrop->SetBinContent(xBinNew, yBinNew, content);
                }
            }
        }
    }
    return hCrop;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <data>" << argv[1] << " <partition>" << argv[2] << " <cell>" << argv[3] << std::endl;
        // data = {0: muon, 1: nue, 2: data}
        // cell for MC is brick
        return 1;
    }
    int data = std::atoi(argv[1]);
    int partition = std::atoi(argv[2]);
    int cell = std::atoi(argv[3]);

    setRange(data, &path, cell, &xMin, &xMax, &yMin, &yMax, &xBin, &yBin, &nPlates);
    TH2::AddDirectory(false);
    TH2F* h2, *hCrop;
    TH3F* h3;
    TH2F* hComb = new TH2F("XYseg", "XYseg", xBin, xMin, xMax, yBin, yMin, yMax);
    TFile* f, *ff[9];

    if (data==0) {
        TString fileName = TString::Format("%s/hist_XYP_muon.root", path.Data());
        f = TFile::Open(fileName);
        h3 = (TH3F*)(f->Get("XYPseg"));
    }
    else if (data==1) {
        TString fileName = TString::Format("%s/hist_XYP_nue.root", path.Data());
        f = TFile::Open(fileName);
    }
    else if (data==2) {
        openFiles(cell, &ff[0]);
    }

    int combination = 0;
    int combStart = partition * 1000;
    int combEnd = combStart + 1000;

    TStopwatch stopWatch;
    stopWatch.Start();
    for (double shiftTX = -shiftRange; shiftTX <= shiftRange; shiftTX += shiftStep) {
        for (double shiftTY = -shiftRange; shiftTY <= shiftRange; shiftTY += shiftStep) {
            stopWatch.Continue();
            combination++;
            if (combination < combStart || combination >= combEnd) continue;

            TString outputPath = TString::Format("%s/../shift/%i", path.Data(), cell);
            if (!std::filesystem::exists(outputPath.Data())) std::filesystem::create_directory(outputPath.Data());
            TString outputName = TString::Format("%s/histo_shifts_%i.root", outputPath.Data(), combination);
	        TFile *outputFile = new TFile(outputName, "RECREATE");
            std::cout << "Combination " << combination << std::endl;

            hComb->Reset();
            for (int layer = 0; layer < nPlates; ++layer) {
                int plate = layer + 1;
                double shiftX = shiftTX / 1000.0 * stepZ * layer;
                double shiftY = shiftTY / 1000.0 * stepZ * layer;

                TString histName = TString::Format("XYseg_%d", plate);
                TString histTitle = TString::Format("XYPseg_%d", plate);

                if (data == 0) {
                    h3->GetZaxis()->SetRange(layer,layer);
                    h2 = (TH2F*)(h3->Project3D("yx"));
                    h2->SetNameTitle(histName.Data(), histTitle);
                }
                else if (data==1) {
                    h2 = (TH2F*)(f->Get(histName.Data()));
                }
                else if (data==2) {
                    h2 = matrixCells(&ff[0], plate);
                }
                hCrop = cropHist(h2, shiftX, shiftY);
                delete h2;
                hComb->Add(hCrop);
                outputFile->cd();
                hCrop->Write();
                delete hCrop;
                // std::cout << plate << std::endl;
            }
            outputFile->cd();
            hComb->Write();
            outputFile->Write();
            outputFile->Close();
            delete outputFile;
            std::cout << "Time: " << stopWatch.RealTime() << std::endl;
            printMemoryInfo();
        }
    }
    delete h3;
    if (data!=2) f->Close();
    else {
        for (int i = 0; i < 9; i++) {
            ff[i]->Close();
        }
    }
    return 0;
}