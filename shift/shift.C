#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TKey.h"
#include "TROOT.h"
#include "TList.h"
#include "TSystem.h"
#include <iostream>
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

// Parameters
const int binSize    = 50;  //um
const int shiftRange = 50; //mrad
const int shiftStep  = 2;   //mrad
const int stepZ      = 1000 + 350; //um
int xMin, xMax, yMin, yMax, xBin, yBin;
TString path;
// TString path = "/Users/fabioali/cernbox"; //for local test

void setRange(int data, TString* path, int fragment, int* xMin, int* xMax, int* yMin, int* yMax, int* xBin, int* yBin) {
    if (data==0) {
        *path = TString::Format("/eos/user/f/falicant/Simulations_sndlhc/muon1E5_simsndlhc/b%05i", fragment);
        *xMin = 284000;
        *xMax = 304000;
        *yMin = 79000;
        *yMax = 99000;
    }
    else if (data==1) {
        *path = TString::Format("/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b%05i", fragment);
        *xMin = 0;
        *xMax = 200000;
        *yMin = 0;
        *yMax = 200000;
    }
    else if (data==2) {
        *path = "/Users/fabioali/cernbox";
        // *path = "/eos/user/f/falicant/nue_search/R1B121/gen1/hist";
        const int xLow = fragment / 18;
        const int yLow = fragment % 18;
        *xMin = xLow*10000 + 3000;
        *xMax = xLow*10000 + 17000;
        *yMin = yLow*10000 + 3000;
        *yMax = yLow*10000 + 17000;
    }
    *xBin = int((*xMax - *xMin) / binSize);
    *yBin = int((*yMax - *yMin) / binSize);
}

void openFiles(int fragment, TFile* f[9]) {
    int idx = 0;
    for (int yCell = fragment-18; yCell <= fragment+18; yCell +=18) {
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
        // std::cout << f[i] << std::endl;
        h2 = (TH2F*)f[i]->Get(histName);
        list->Add(h2);
    }
    hm->Merge(list);
    delete h2;
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
        std::cerr << "Usage: " << argv[0] << " <data>" << argv[1] << " <partition>" << argv[2] << " <fragment>" << std::endl;
        // data = {0: muon, 1: nue, 2: data}
        // fragment for MC is brick
        return 1;
    }
    int data = std::atoi(argv[1]);
    int partition = std::atoi(argv[2]);
    int fragment = std::atoi(argv[3]);

    setRange(data, &path, fragment, &xMin, &xMax, &yMin, &yMax, &xBin, &yBin);
    
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
        openFiles(fragment, &ff[0]);
    }

    int combination = 0;
    int combStart = partition * 100;
    int combEnd = combStart + 5;
    
    for (double shiftTX = -shiftRange; shiftTX <= shiftRange; shiftTX += shiftStep) {
        for (double shiftTY = -shiftRange; shiftTY <= shiftRange; shiftTY += shiftStep) {
            combination++;
            if (combination < combStart || combination >= combEnd) continue;

            TString outputName = TString::Format("%s/shift/histo_shifts_%i_%i.root", path.Data(), fragment, combination);
	        TFile *outputFile = new TFile(outputName, "RECREATE");
            std::cout << "Combination " << combination << std::endl;

            hComb->Reset();
            for (int layer = 0; layer < 57; ++layer) {
                int plate = layer + 1;
                double shiftX = shiftTX / 1000.0 * stepZ * layer;
                double shiftY = shiftTY / 1000.0 * stepZ * layer;

                TString histName = TString::Format("XYseg_%d", plate);
                TString histTitle = TString::Format("XYPseg_%d", plate);

                if (data == 0) {
                    h3->GetZaxis()->SetRange(plate,plate);
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
                std::cout << plate << std::endl;
            }
            outputFile->cd();
            hComb->Write();
            outputFile->Close();
            delete outputFile;
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


// #include <stdio.h>

// // Function to determine the group of the brick
// int getBrickGroup(int brick) {
//     int bricks1[] = {11, 21, 31, 41, 51};
//     int bricks2[] = {12, 22, 32, 42, 52};
//     int bricks3[] = {13, 23, 33, 43, 53};
//     int bricks4[] = {14, 24, 34, 44, 54};
    
//     // Check if the brick belongs to any of the groups
//     for (int i = 0; i < 5; i++) {
//         if (brick == bricks1[i]) return 1;
//         if (brick == bricks2[i]) return 2;
//         if (brick == bricks3[i]) return 3;
//         if (brick == bricks4[i]) return 4;
//     }
    
//     return 0;  // Return 0 if brick is not found
// }

// int main() {
//     int brick = 31;  // Example brick
//     int offsetx = 0;
//     int offsety = 0;
    
//     // Get the group of the brick
//     int group = getBrickGroup(brick);

//     // Switch based on the group number
//     switch (group) {
//         case 1:
//             offsetx = 200000;
//             offsety = 0;
//             break;
//         case 2:
//             offsetx = 0;
//             offsety = 0;
//             break;
//         case 3:
//             offsetx = 200000;
//             offsety = 200000;
//             break;
//         case 4:
//             offsetx = 0;
//             offsety = 200000;
//             break;
//         default:
//             printf("Brick not found in any group\n");
//             return 1;  // Exit with error code
//     }
    
//     printf("OffsetX: %d, OffsetY: %d\n", offsetx, offsety);
//     return 0;
// }