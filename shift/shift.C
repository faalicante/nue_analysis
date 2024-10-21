#include "TFile.h"
#include "TH2F.h"
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
    while ((key = (TKey*)next())) {
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

// const int xMin = 0;       //add switch case (at the end of the code)
// const int xMax = 200000;
// const int yMin = 0;
// const int yMax = 200000;

// Muons
const int xMin = 289000;
const int xMax = 300000;
const int yMin = 84000;
const int yMax = 95000;

const int xBin = int((xMax - xMin) / binSize);
const int yBin = int((yMax - yMin) / binSize);

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
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <partition>" << std::endl;
        return 1;
    }
    int partition = std::atoi(argv[1]);

    // const char* path = "/Users/fabioali/cernbox";
    // const char* path = "/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b000022";
    const char* path = "/eos/user/f/falicant/Simulations_sndlhc/muon1E5_simsndlhc/b000021";
    // TString file = TString::Format("%s/hist_XYP_nue.root", path);
    TString file = TString::Format("%s/hist_XYP_muon.root", path);
    // std::map<std::string, TH2F*> h = loadHists(file.Data());

    int combination = 0;
    int combStart = partition * 100;
    int combEnd = combStart + 100;

    for (double shiftTX = -shiftRange; shiftTX <= shiftRange; shiftTX += shiftStep) {
        for (double shiftTY = -shiftRange; shiftTY <= shiftRange; shiftTY += shiftStep) {
            combination++;
            if (combination < combStart || combination >= combEnd) continue;

            TString outputFileName = TString::Format("%s/shift_full/histo_shifts_%d.root", path, combination);
	        TFile *outputFile = new TFile(outputFileName, "RECREATE");
            std::cout << "Combination " << combination << std::endl;

            TH2F* hComb = new TH2F("XYseg", "XYseg", xBin, xMin, xMax, yBin, yMin, yMax);
            
            for (int layer = 0; layer < 60; ++layer) {
                int plate = layer + 1;
                double shiftX = shiftTX / 1000.0 * stepZ * layer;
                double shiftY = shiftTY / 1000.0 * stepZ * layer;
                
                // TString histName = TString::Format("XYseg_%d", plate);
                TH3F *h3 = (TH3F *)(file->Get("XYPseg"));
                h3->GetZaxis()->SetRange(plate,plate);
                TH2F *h = (TH2F*)(h3->Project3D("yx"));
                TH2F* hCrop = cropHist(h, shiftX, shiftY);
                // TH2F* hCrop = cropHist(h[histName.Data()], shiftX, shiftY);
                hCrop->Write();
                hComb->Add(hCrop);
                delete hCrop;
            }
            hComb->Write();
            delete hComb;
            outputFile->Close();
	    delete outputFile;
        }
    }

    // Clean up
    // for (auto& pair : h) {
    //     delete pair.second;
    // }

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