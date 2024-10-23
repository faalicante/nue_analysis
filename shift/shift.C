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
        std::cout << key->GetName() << std::endl;
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
// int xMin = 289000;
// int xMax = 300000;
// int yMin = 84000;
// int yMax = 95000;

// int xBin = int((xMax - xMin) / binSize);
// int yBin = int((yMax - yMin) / binSize);

const char* path = "/eos/user/f/falicant/nue_search/R1B121/gen1/hist";
// const char* path = "/Users/fabioali/cernbox";


void setRange(int fragment, int* xLow, int* yLow) {
    TString histFile = TString::Format("%s/hist_XYP_b121_%i.root", path, fragment);
    TFile* f = TFile::Open(histFile);
    TH2F* h2 = (TH2F*)f->Get("XYseg_1");
    *xLow = h2->GetXaxis()->GetXmin();
    *yLow = h2->GetYaxis()->GetXmin();
    f->Close();
}

TH2F* matrixCells(int fragment, int plate) {
    TList *list = new TList;
    // TFile* f;
    // TH2F* h2;
    TString histName = TString::Format("XYseg_%d", plate);
    for (int yCell = fragment-18; yCell <= fragment+18; yCell +=18) {
        for (int xCell = yCell-1; xCell <= yCell+1; xCell++) {
            TString histFile = TString::Format("%s/hist_XYP_b121_%i.root", path, xCell);
            // f = TFile::Open(histFile);
            // TH2F* h2 = (TH2F*)f->Get(histName.Data());
            // h2->SetDirectory(gROOT);
            // h2->SetName(histName);
            
            std::map<std::string, TH2F*> h = loadHists(histFile.Data());
            list->Add(h[histName.Data()]);
            // list->Add(h2);
            // delete h2;
        }
    }
    
    TH2F *hm = (TH2F*)(((TH2F*)list->At(0))->Clone("hm"));
    hm->Reset();
    hm->Merge(list);
    delete list;
    // f->Close();

    return hm;
}

TH2F* cropHist(TH2F* h2, double shiftX, double shiftY, int xBin, int xMin, int xMax, int yBin, int yMin, int yMax) {
    std::cout << h2->GetTitle() << std::endl;
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
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <partition>" << std::endl;
        return 1;
    }
    int fragment = std::atoi(argv[1]);
    int partition = std::atoi(argv[2]);

    int xLow = 0;
    int yLow = 0;


    TString histFile = TString::Format("%s/hist_XYP_b121_%i.root", path, fragment);
    // std::map<std::string, TH2F*> h = loadHists(file.Data());
    setRange(fragment, &xLow, &yLow);
    // setRange(h["XYseg_1"], &xLow, &yLow);
    int xMax = xLow + 11000;
    int yMax = yLow + 11000;
    int xMin = xLow - 1000;
    int yMin = yLow - 1000;
    int xBin = int((xMax - xMin) / binSize);
    int yBin = int((yMax - yMin) / binSize);
    
    // const char* path = "/Users/fabioali/cernbox";
    // const char* path = "/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b000022";
    // const char* path = "/eos/user/f/falicant/Simulations_sndlhc/muon1E5_simsndlhc/b000021";
    // TString file = TString::Format("%s/hist_XYP_nue.root", path);
    // TFile *inputFile = TFile::Open(TString::Format("%s/hist_XYP_muon.root", path));

    int combination = 0;
    int combStart = partition * 100;
    int combEnd = combStart + 100;

    for (double shiftTX = -shiftRange; shiftTX <= shiftRange; shiftTX += shiftStep) {
        for (double shiftTY = -shiftRange; shiftTY <= shiftRange; shiftTY += shiftStep) {
            combination++;
            if (combination < combStart || combination >= combEnd) continue;

            TString outputFileName = TString::Format("%s/shift_full/histo_shifts_%i_%i.root", path, fragment, combination);
	        TFile *outputFile = new TFile(outputFileName, "RECREATE");
            std::cout << "Combination " << combination << std::endl;

            // TList *list = new TList;
            TH2F* hComb = new TH2F("XYseg", "XYseg", xBin, xMin, xMax, yBin, yMin, yMax);
            
            // TH2F* hm;
            // TH2F* hCrop;
            for (int layer = 0; layer < 57; ++layer) {
                int plate = layer + 1;
                double shiftX = shiftTX / 1000.0 * stepZ * layer;
                double shiftY = shiftTY / 1000.0 * stepZ * layer;
                
                std::cout << plate << std::endl;
                // TH3F *h3 = (TH3F*)(inputFile->Get("XYPseg"));
                // h3->GetZaxis()->SetRange(plate,plate);
                // TH2F *h = (TH2F*)(h3->Project3D("yx"));
                // h->SetTitle(Form("Plate %d",plate));
                // TH2F* hCrop = cropHist(h, shiftX, shiftY);
                // TString histName = TString::Format("XYseg_%d", plate);
                TH2F* hm = matrixCells(fragment, plate);
                // hm->SetName(histName.Data());
                TH2F* hCrop = cropHist(hm, shiftX, shiftY, xBin, xMin, xMax, yBin, yMin, yMax);
                // TH2F* hCrop = cropHist(h[histName.Data()], shiftX, shiftY);
                outputFile->cd();
                hCrop->Write();
                // list->Add(hCrop);
                hComb->Add(hCrop);
                delete hm;
                delete hCrop;
            }
            // TH2F *hComb = (TH2F*)(((TH2F*)list->At(0))->Clone("hComb"));
            // hComb->Reset();
            // hComb->Merge(list);
            // hComb->SetNameTitle("XYseg","XYseg");
            // delete list;
            outputFile->cd();
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

// make case for sig bkg data
