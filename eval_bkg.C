#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2F.h>
#include <TString.h>

void eval_bkg(int cell, float x0, float y0) {
    TString path = TString("/eos/experiment/sndlhc/users/falicant/RUN1/b121/hist");
    TString fname = TString::Format("%s/hist_XYP_b121_%i.root", path.Data(), cell);
    TFile* f = TFile::Open(fname);
    TH2F *h2 = (TH2F *)(f->Get("XYseg"));
    TFile *file = new TFile("bkg_b121.root", "UPDATE");
    TTree *tree = (TTree*)file->Get("bkg");
    
    float x;
    float y;
    float b;

    if (!tree) {
        tree = new TTree("bkg", "Segment background");
        tree->Branch("x", &x);
        tree->Branch("y", &y);
        tree->Branch("b", &b);
    }
    else {
        tree->SetBranchAddress("x", &x);
        tree->SetBranchAddress("y", &y);
        tree->SetBranchAddress("b", &b);
    }

    int r0 = (int)(500./50);
    x0 = x0 * 1000;
    y0 = y0 * 1000;
    int xbin = h2->GetXaxis()->FindBin(x0);
    int ybin = h2->GetYaxis()->FindBin(y0);

    h2->Smooth();
    for (int ix = xbin-r0; ix < xbin+r0; ++ix) {
        for (int iy = ybin-r0; iy < ybin+r0; ++iy) {
            float x1 = h2->GetXaxis()->GetBinLowEdge(ix);
            float y1 = h2->GetYaxis()->GetBinLowEdge(iy);
            int b1 = h2->GetBinContent(ix, iy);
            printf("x: %f, y: %f, b: %i\n",x1,y1,b1);
            x = x1;
            y = y1;
            b = b1;
            tree->Fill();
        }
    }
    
    file->cd();
    tree->Write("", TObject::kOverwrite);
    file->Close();
}