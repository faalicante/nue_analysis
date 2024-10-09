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

// Function to crop histogram
TH2F* cropHist(TH2F* h2, double shiftX, double shiftY) {
    const int xBin = (82500 - 71500) / 50;
    const int yBin = (96500 - 85500) / 50;
    TH2F* hCrop = new TH2F(h2->GetTitle(), h2->GetTitle(), xBin, 71500, 82500, yBin, 85500, 96500);
    for (int xBin = 1; xBin <= h2->GetNbinsX(); ++xBin) {
        double xCenter = h2->GetXaxis()->GetBinCenter(xBin) + shiftX;
        if (xCenter <= 82500 && xCenter >= 71500) {
            for (int yBin = 1; yBin <= h2->GetNbinsY(); ++yBin) {
                double yCenter = h2->GetYaxis()->GetBinCenter(yBin) + shiftY;
                if (yCenter <= 96500 && yCenter >= 85500) {
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

    const char* path = "/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b000022";
    TString file = TString::Format("%s/hist_XYP_nue.root", path);
    std::map<std::string, TH2F*> h = loadHists(file.Data());

    int combination = 0;
    int combStart = partition * 100;
    int combEnd = combStart + 100;

    for (double shiftTX = -30; shiftTX <= 30; shiftTX += 2) {
        for (double shiftTY = -30; shiftTY <= 30; shiftTY += 2) {
            combination++;
            if (combination < combStart || combination >= combEnd) continue;

            TString outputFileName = TString::Format("%s/histo_shifts_%d.root", path, combination);
	    //For testing
//            TString outputFileName = Form("histo_shifts_%d.root",combination);
//	    TFile *outputFile = new TFile("/eos/user/m/maclimes/SND/SND_nue_data/simulation/JT_Fabio/test/"+outputFileName,"RECREATE");	
	    TFile *outputFile = new TFile(outputFileName, "RECREATE");
            std::cout << "Combination " << combination << std::endl;

            TH2F* hComb = new TH2F("XYseg", "XYseg", (82500 - 71500) / 50, 71500, 82500, (96500 - 85500) / 50, 85500, 96500);
            
            for (int layer = 0; layer < 60; ++layer) {
                int plate = layer + 1;
                double shiftX = shiftTX / 1000.0 * 1350 * layer;
                double shiftY = shiftTY / 1000.0 * 1350 * layer;
                
                TString histName = TString::Format("XYseg_%d", plate);
                TH2F* hCrop = cropHist(h[histName.Data()], shiftX, shiftY);
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
    for (auto& pair : h) {
        delete pair.second;
    }

    return 0;
}
