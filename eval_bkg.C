#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2F.h>
#include <TString.h>

int r0 = (int)(250./50);
float x;
float y;
int p;
int b;

TH2F* projectHist(TFile* f, int plate) {
    TH3F* h3 = (TH3F*)(f->Get("XYPseg"));
    h3->GetZaxis()->SetRange(plate+1,plate+1);
    TH2F* h2 = (TH2F*)(h3->Project3D("yx"));
    return h2;
}

void fillTree(TH2F* h2, TTree* tree, float x0, float y0, int plate) {
    int xbin = h2->GetXaxis()->FindBin(x0);
    int ybin = h2->GetYaxis()->FindBin(y0);
    h2->Smooth();
    for (int ix = xbin-r0; ix < xbin+r0; ++ix) {
        for (int iy = ybin-r0; iy < ybin+r0; ++iy) {
            float x1 = h2->GetXaxis()->GetBinLowEdge(ix);
            float y1 = h2->GetYaxis()->GetBinLowEdge(iy);
            int b1 = h2->GetBinContent(ix, iy);
            printf("x: %f, y: %f, p: %i, b: %i\n",x1,y1,plate,b1);
            x = x1;
            y = y1;
            p = plate;
            b = b1;
            tree->Fill();
        }
    }
}

void eval_bkg(int cell, float x0, float y0) {
    const int xCell = cell % 18;
    const int yCell = cell / 18;
    // TString path = TString("/Users/fabioali/cernbox/test_shift/trk");
    TString path = TString("/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b000121/cells");
    TString fname = TString::Format("%s/cell_%i0_%i0_1x1cm/b000021/b000021.0.0.0.trk.root", path.Data(), xCell+1, yCell+1);
    TFile* f = TFile::Open(fname);
    TFile *file = new TFile("bkg_b121.root", "UPDATE");
    TTree *tree = (TTree*)file->Get("bkg");
    TH2F *h2;

    if (!tree) {
        tree = new TTree("bkg", "Segment background");
        tree->Branch("x", &x);
        tree->Branch("y", &y);
        tree->Branch("p", &p);
        tree->Branch("b", &b);
    }
    else {
        tree->SetBranchAddress("x", &x);
        tree->SetBranchAddress("y", &y);
        tree->SetBranchAddress("p", &p);
        tree->SetBranchAddress("b", &b);
    }

    x0 = x0 * 1000;
    y0 = y0 * 1000;

    h2 = (TH2F *)(f->Get("XYseg"));
    fillTree(h2, tree, x0, y0, 0);
    delete h2;
    for (p=1; p<=57; p++) {
        h2 = projectHist(f, p);
        fillTree(h2, tree, x0, y0, p);
        delete h2;
    }

    file->cd();
    tree->Write("", TObject::kOverwrite);
    file->Close();
}