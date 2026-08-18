#include "TFile.h"
#include "TTree.h"
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
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <limits>
#include <string>

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
const int radius     = 500;  // (um)
const int ntag = 250;
const int nPlates = 57;
const int cropSize = 20;
const int hardZScaleSize = 10;
const int cropVoxels = nPlates * cropSize * cropSize;
const int dz = 1350;
const Int_t poissonSample = 0;
const Int_t hardNegativeSample = 1;
const Int_t signalSample = 2;
int xMin, xMax, yMin, yMax, xBins, yBins, xLow, yLow;
int range;
float bkg = 0.0F;
double threshold = 0.0;
TString path;
TString opath;
TString ppath;
TString histName;

class SampleTreeWriter {
public:
    explicit SampleTreeWriter(const TString& outputPath) {
        outputFile_ = TFile::Open(outputPath, "RECREATE");
        tree_ = new TTree("samples", "Raw source regions");
        tree_->SetDirectory(outputFile_);
        const TString countsLeaf = TString::Format(
            "counts[%d][%d][%d]/I", nPlates, cropSize, cropSize);
        tree_->Branch("counts", counts_, countsLeaf.Data());
        tree_->Branch("background_mu", &backgroundMu_, "background_mu/F");
        tree_->Branch("presence", &presence_, "presence/I");
        tree_->Branch("sample_type", &sampleType_, "sample_type/I");
        tree_->Branch("slope_x", &slopeX_, "slope_x/F");
        tree_->Branch("slope_y", &slopeY_, "slope_y/F");
        tree_->Branch("signal_event_id", &signalEventId_, "signal_event_id/I");
        tree_->Branch("cell_id", &cellId_, "cell_id/I");
        tree_->Branch("tag_cell_id", &tagCellId_, "tag_cell_id/I");
    }

    void fill(TH2** rawLayers, float x0, float y0, int presence, int sampleType,
              float slopeX, float slopeY, int signalEventId,
              int cellId, int tagCellId, float bkg) {
        std::fill_n(counts_, cropVoxels, 0);
        presence_ = presence;
        sampleType_ = sampleType;
        slopeX_ = slopeX;
        slopeY_ = slopeY;
        signalEventId_ = signalEventId;
        cellId_ = cellId;
        tagCellId_ = tagCellId;
        backgroundMu_ = bkg/nPlates;

        for (int z = 0; z < nPlates; ++z) {
            TH2* layer = rawLayers[z];

            const int centerX = layer->GetXaxis()->FindFixBin(x0);
            const int centerY = layer->GetYaxis()->FindFixBin(y0);
            const int firstX = centerX - cropSize / 2;
            const int firstY = centerY - cropSize / 2;

            for (int y = 0; y < cropSize; ++y) {
                const int sourceY = firstY + y;
                for (int x = 0; x < cropSize; ++x) {
                    const int sourceX = firstX + x;
                    int value = 0;
                    if (sourceX > 0 && sourceX <= layer->GetNbinsX() &&
                        sourceY > 0 && sourceY <= layer->GetNbinsY()) {
                        value = layer->GetBinContent(sourceX, sourceY);
                    }
                    const int index = (z * cropSize * cropSize) + (y * cropSize) + x;
                    counts_[index] = value;
                }
            }
        }
        tree_->Fill();
    }

    void write() {
        outputFile_->cd();
        tree_->Write();
        outputFile_->Close();
        delete outputFile_;
        outputFile_ = nullptr;
        tree_ = nullptr;
    }

private:
    TFile* outputFile_ = nullptr;
    TTree* tree_ = nullptr;
    Int_t counts_[cropVoxels] = {};
    Float_t backgroundMu_ = 0.0F;
    Int_t presence_ = 0;
    Int_t sampleType_ = poissonSample;
    Float_t slopeX_ = 0.0F;
    Float_t slopeY_ = 0.0F;
    Int_t signalEventId_ = -1;
    Int_t cellId_ = -1;
    Int_t tagCellId_ = -1;
};

void getPath(int data, TString* path, TString* opath, TString *ppath, int cell, int* xLow, int* yLow, int* range) {
    if (data == 0 || data == 3) { // Muon simulation //data==3 is for the "nothing" sample
        // *path = "/Users/fabioali/cernbox/shift/muon";
        // *opath = *path;
        *path = TString::Format("/eos/experiment/sndlhc/users/dancc/FEDRA/muon_Euniform_RUN1_FLUKA25/b%06i/cell_reco", brick);
        *opath = "/eos/experiment/sndlhc/users/falicant/RNN2/none";
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
        *opath = "/eos/experiment/sndlhc/users/falicant/RNN2/signal";
        *ppath = TString::Format("%s/%i", opath->Data(), cell);
        *range = 500;
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
        *range = 500;
    }
}

TH2F* setRanges1(int data, int cell, TFile** f, int* xMin, int* xMax, int* yMin, int* yMax, int* xBins, int* yBins) {
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
    h2->GetXaxis()->SetRangeUser(*xMin, *xMax);
    h2->GetYaxis()->SetRangeUser(*yMin, *yMax);
    return h2;
}

void setRanges2(TH2 **hm, float x0, float y0) {
    for (int i = 0; i < nPlates; ++i) {
        hm[i]->GetXaxis()->SetRangeUser(x0-range, x0+range);
        hm[i]->GetYaxis()->SetRangeUser(y0-range, y0+range);
    }
    // std::cout << "x range = [" << x0-range << ", " << x0+range << "], y range = [" << y0-range << ", " << y0+range << "]" << std::endl;
}

TH3F* loadH3(TFile *f) {
    TH3F *h3 = nullptr;
    if (f) {
        f->GetObject("XYPseg", h3);
        if (h3) h3->SetDirectory(0);
    }
    return h3;
}

TH1F* drawSpectrum(TH2F *h2, const char* name) {
    int nBinsX = h2->GetNbinsX();
    int nBinsY = h2->GetNbinsY();
    TH1F* hSpec = new TH1F(name, "Spectrum;rankbin", 200, 0, 1000);
    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            int content = h2->GetBinContent(i, j);
            if (content > 0) hSpec->Fill(content);
        }
    }
    return hSpec;
}

void poisBkg(TH1F* rawSpectrum, TH1F* smoothedSpectrum,
             float *bkg, double *threshold) {
    const double rawMpv = rawSpectrum->GetBinLowEdge(rawSpectrum->GetMaximumBin());
    const double smoothedMpv =
        smoothedSpectrum->GetBinLowEdge(smoothedSpectrum->GetMaximumBin());
    *bkg = static_cast<float>(rawMpv);
    *threshold = smoothedMpv + 3*std::sqrt(smoothedMpv);
    std::cout << "Poisson background (raw MPV): " << *bkg << std::endl;
    std::cout << "Poisson threshold (smoothed): " << *threshold << std::endl;
}

void openFiles(int data, int cell, const TString& inputPath, TFile** f, TH3F** H3cell) {
    TString fileName;
    if (!inputPath.IsNull()) fileName = inputPath;
    else if (data==1) fileName = TString::Format("%s/b%06i.0.0.%i.trk.root", path.Data(), brick, cell+1);
    else if (data==0 || data==3) fileName = TString::Format("%s/cell_%i0_%i0/b%06i/b%06i.0.%i.%i.trk.root", path.Data(), xLow, yLow, brick, brick, xLow, yLow);
    // std::cout << fileName << std::endl;
    *f = TFile::Open(fileName);
    TH2F* H2cell = setRanges1(data, cell, f, &xMin, &xMax, &yMin, &yMax, &xBins, &yBins);
    *H3cell = loadH3(*f);
    TH1F* rawSpectrum = drawSpectrum(H2cell, "raw_background_spectrum");
    TH2F* smoothedH2cell = (TH2F*)H2cell->Clone("XYseg_background_smoothed");
    smoothedH2cell->SetDirectory(nullptr);
    smoothedH2cell->Smooth();
    TH1F* smoothedSpectrum =
        drawSpectrum(smoothedH2cell, "smoothed_background_spectrum");
    poisBkg(rawSpectrum, smoothedSpectrum, &bkg, &threshold);
    delete rawSpectrum;
    delete smoothedSpectrum;
    delete smoothedH2cell;
}

TH2* projectHist(TH3F* h3, int plate) {
    h3->GetZaxis()->SetRange(plate+1,plate+1);
    TH2* h2 = (TH2*)(h3->Project3D("yx"));
    return h2;
}

TH2* stackHist(TH2 **rawLayers, TH2 **smoothedLayers, TH3F *H3cell) {
    TH2* hComb = nullptr;
    for (int layer = 0; layer < nPlates; ++layer) {
        int plate = layer + 1;
        rawLayers[layer] = projectHist(H3cell, plate);
        if (hComb == nullptr) {
            hComb = (TH2*)rawLayers[layer]->Clone("XYseg_candidates");
            hComb->Reset();
        }
        hComb->Add(rawLayers[layer]);
        smoothedLayers[layer] = (TH2*)rawLayers[layer]->Clone(
            TString::Format("XYseg_smoothed_%d", plate));
        smoothedLayers[layer]->Smooth();
    }
    hComb->Smooth();
    return hComb;
}

void getMax(int data, TH2 &h2, TObjArray &peaks, double threshold) {
    const double rankbin = h2.GetMaximum();
    Int_t MaxBin = h2.GetMaximumBin();
    Int_t ix,iy,iz;
    h2.GetBinXYZ(MaxBin, ix, iy, iz);
    float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
    float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
    TEllipse  *el = new TEllipse(x,y,radius,radius);
    int r0 = (int)round((double)radius/binSize);
    for(int iix = ix-r0; iix<=ix+r0; iix++) {
        for(int iiy = iy-r0; iiy<=iy+r0; iiy++) {
            double dx = iix - ix;
            double dy = iiy - iy;
            double distance = (dx*dx + dy*dy)/(r0*r0);
            if (distance <= 1) h2.SetBinContent(iix,iiy,0);
        }
    }
    if ( (data == 0 && rankbin > threshold) || (data == 3 && rankbin < threshold) ) {
        el->SetFillStyle(0);
        peaks.Add(el);
    }
    else {
        delete el;
    }
}

void get_peaks(int data,TH2 &h2, TObjArray &peaks, int npmax, double threshold) {
    TH2 *h2new = (TH2*)h2.Clone("get_peaks");
    for(int i=0; i<npmax; i++){
        getMax(data,*h2new, peaks, threshold);
    }
    delete h2new;
}

double findColScale(TH2 **hm, float x0, float y0, int windowSize) {
    double maxValue = 0;
    for (int i = 0; i < nPlates; ++i) {
        if (hm[i] == nullptr) continue;
        const int centerX = hm[i]->GetXaxis()->FindFixBin(x0);
        const int centerY = hm[i]->GetYaxis()->FindFixBin(y0);
        const int firstX = centerX - windowSize / 2;
        const int firstY = centerY - windowSize / 2;
        for (int y = 0; y < windowSize; ++y) {
            for (int x = 0; x < windowSize; ++x) {
                const double value = hm[i]->GetBinContent(firstX + x, firstY + y);
                if (value > maxValue) maxValue = value;
            }
        }
    }
    // std::cout << "Max value across all plates: " << maxValue << std::endl;
    return maxValue;
}

void printBW (TH2 **hm, int cell, double zScale, int tag) {
    if (!std::filesystem::exists(ppath.Data())) std::filesystem::create_directory(ppath.Data());
    if (!std::filesystem::exists(TString::Format("%s/%i", ppath.Data(), tag).Data())&&tag>0) {
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

void setStyle() {
    gErrorIgnoreLevel = kWarning;
    gROOT->SetBatch(!print);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(52);
    gStyle->SetTitleSize(0, "XYZ");
    gStyle->SetLabelSize(0, "XYZ");
    gStyle->SetTickLength(0, "XYZ");
    gStyle->SetFrameLineWidth(0);
    gStyle->SetPadLeftMargin(0);
    gStyle->SetPadRightMargin(0);
    gStyle->SetPadTopMargin(0);
    gStyle->SetPadBottomMargin(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetCanvasColor(1);
    gStyle->SetPadColor(1);
    gROOT->ForceStyle();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <data> <cell> [--input file.root] [--output samples.root]"
                  << " [--x0 X --y0 Y --tx TX --ty TY --p0 PLATE] [--images]" << std::endl;
        // data = {0: muon, 1: nue, 2: data}
        // cell for nue is event
        return 1;
    }
    int data = std::atoi(argv[1]);
    int cell = std::atoi(argv[2]);
    
    // for neutrino only
    double xn = std::numeric_limits<double>::quiet_NaN();
    double yn = std::numeric_limits<double>::quiet_NaN();
    double txn = std::numeric_limits<double>::quiet_NaN();
    double tyn = std::numeric_limits<double>::quiet_NaN();
    int pn = -1;
    TString inputPath;
    TString outputPath = "samples.root";
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--x0" && i + 1 < argc) {
            xn = std::stod(argv[++i]);
        }
        else if (arg == "--y0" && i + 1 < argc) {
            yn = std::stod(argv[++i]);
        }
        else if (arg == "--tx" && i + 1 < argc) {
            txn = std::stod(argv[++i]);
        }
        else if (arg == "--ty" && i + 1 < argc) {
            tyn = std::stod(argv[++i]);
        }
        else if (arg == "--p0" && i + 1 < argc) {
            pn = std::atoi(argv[++i]);
        }
        else if (arg == "--input" && i + 1 < argc) {
            inputPath = argv[++i];
        }
        else if (arg == "--output" && i + 1 < argc) {
            outputPath = argv[++i];
        }
        else if (arg == "--images") {
            print = true;
        }
        else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            return 1;
        }
    }
    if (data == 1 && (!std::isfinite(xn) || !std::isfinite(yn) ||
                      !std::isfinite(txn) || !std::isfinite(tyn) || pn < 0)) {
        std::cerr << "data=1 requires --x0, --y0, --tx, --ty and --p0" << std::endl;
        return 1;
    }
    
    TStopwatch stopWatch;
    stopWatch.Start();

    setStyle();
    
    getPath(data, &path, &opath, &ppath, cell, &xLow, &yLow, &range);
    
    TFile *f = nullptr;
    TH3F *H3cell = nullptr;

    openFiles(data, cell, inputPath, &f, &H3cell);

    TH2 *rawLayers[nPlates] = {};
    TH2 *smoothedLayers[nPlates] = {};
    TH2::AddDirectory(false);
    TH2 *hComb = stackHist(&rawLayers[0], &smoothedLayers[0], H3cell);
    TObjArray peaks;
    peaks.SetOwner(kTRUE);
    get_peaks(data,*hComb,peaks,ntag,threshold);
    SampleTreeWriter sampleWriter(outputPath);
    int outputSamples = 0;
    
    if (data == 0 || data == 3) {
        int np = peaks.GetEntries();
        for(int i=0; i<np; i++) {
            TEllipse *el = ((TEllipse*)(peaks.At(i)));
            float x0 = el->GetX1();
            float y0 = el->GetY1();
            if ((xLow+19)*10000 > x0-range || x0+range > (xLow+20)*10000 ||
                (yLow-1+0.45)*10000 > y0-range || (y0+range > (yLow+0.45)*10000))
                continue;
            setRanges2(&smoothedLayers[0], x0, y0);
            const int zScaleWindow = data == 3 ? cropSize : hardZScaleSize;
            double zScale = findColScale(&smoothedLayers[0], x0, y0, zScaleWindow);

            if (data == 0 && zScale > 3*(threshold/nPlates)) {
                // std::cout << "zScale = " << zScale << std::endl;
                    sampleWriter.fill(&rawLayers[0], x0, y0, 1, hardNegativeSample,
                                      0.0, 0.0, -1, cell, i+1, bkg);
                    ++outputSamples;
            }
            else if (data == 3 && zScale < 2*(threshold/nPlates)) {
                    sampleWriter.fill(&rawLayers[0], x0, y0, 0, poissonSample,
                                      0.0, 0.0, -1, cell, i+1, bkg);
                    ++outputSamples;
            }
            if (print) printBW(&smoothedLayers[0], cell, zScale, i+1);
        }
    }
    
    else if (data == 1) {
        float x0 = xn+txn*dz*0.5*(nPlates-pn)/1000.0;
        float y0 = yn+tyn*dz*0.5*(nPlates-pn)/1000.0;
        sampleWriter.fill(&rawLayers[0], x0, y0, 1, signalSample,
                          txn, tyn, cell, -1, -1, bkg);
        ++outputSamples;
        if (print) {
            setRanges2(&smoothedLayers[0], x0, y0);
            double zScale = findColScale(&smoothedLayers[0], x0, y0, cropSize);
            // std::cout << "Event " << cell << ": zScale = " << zScale << std::endl;
            printBW(&smoothedLayers[0], cell, zScale, 0);
        }
    }

    sampleWriter.write();
    std::cout << "Output entries: " << outputSamples << " in " << outputPath << std::endl;

    // if (!std::filesystem::exists(ppath.Data())) std::filesystem::create_directory(ppath.Data());

    delete hComb;
    for(int p=1; p<=nPlates; p++) {
        delete rawLayers[p-1];
        delete smoothedLayers[p-1];
    }
    std::cout << "---------------------" << std::endl;
    
    f->Close();
    delete H3cell;

    std::cout << "Time: " << round(stopWatch.RealTime()) << std::endl;
    printMemoryInfo();

    return 0;
}
