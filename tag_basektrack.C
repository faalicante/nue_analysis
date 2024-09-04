float set_range(TH3F &h3, TH2F &h2, float thr) {
  int xmi = h2.FindFirstBinAbove(thr,1);
  int xma = h2.FindLastBinAbove(thr,1);
  int ymi = h2.FindFirstBinAbove(thr,2);
  int yma = h2.FindLastBinAbove(thr,2);
  h3.GetXaxis()->SetRange(xmi,xma);
  h3.GetYaxis()->SetRange(ymi,yma);
  h2.GetXaxis()->SetRange(xmi,xma);
  h2.GetYaxis()->SetRange(ymi,yma);
  return h2.Integral()/(xma-xmi+1)/(yma-ymi+1);
}

float set_limits(TH3F &h3, TH2F &h2) {
  float meanbin = h2.Integral()/h2.GetNbinsX()/h2.GetNbinsY();
  float newmean = set_range(h3, h2, meanbin/2);
  float lastmean = set_range(h3, h2, newmean/2);
  printf("meanbin = %f    newmean = %f   lastmean = %f\n",meanbin, newmean, lastmean);
  return lastmean;
}

void eval_bkg(TH2F &h2, TEllipse &el, float *bkg) { //two control regions next to the circle
  float x0 = el.GetX1();
  float y0 = el.GetY1();
  float r1 = el.GetR1();

  float x1 = x0 - r1 * 1.2;
  float x2 = x0 + r1 * 1.2;
  
  int xmi1 = (int)(x1-r1);
  int xma1 = (int)(x1+r1);
  int ymi1 = (int)(y0-r1);
  int yma1 = (int)(y0+r1);

  int xmax1 = h2.GetXaxis()->FindBin(xma1);
  int xmin1 = h2.GetXaxis()->FindBin(xmi1);
  int ymax1 = h2.GetYaxis()->FindBin(yma1);
  int ymin1 = h2.GetYaxis()->FindBin(ymi1);

  int xmi2 = (int)(x2-r1);
  int xma2 = (int)(x2+r1);
  int ymi2 = (int)(y0-r1);
  int yma2 = (int)(y0+r1);

  int xmax2 = h2.GetXaxis()->FindBin(xma2);
  int xmin2 = h2.GetXaxis()->FindBin(xmi2);
  int ymax2 = h2.GetYaxis()->FindBin(yma2);
  int ymin2 = h2.GetYaxis()->FindBin(ymi2);

  float bkg1 = h2.Integral(xmin1,xmax1,ymin1,ymax1)/(xmax1-xmin1+1)/(ymax1-ymin1+1);
  float bkg2 = h2.Integral(xmin2,xmax2,ymin2,ymax2)/(xmax2-xmin2+1)/(ymax2-ymin2+1);
  *bkg = 0.5 * (bkg1 + bkg2);
  // printf("Average background %f\n", *bkg);
}

void eval_bkg(TH2F &h2, float *bkg) { //just one clear area
  float x1 = 292000; //b11 293000 //b21 292000 //b31 292000 //b41 292000 //b51 294000
  float y1 = 92000;  //b11 88000  //b21 91500  //b31 88000  //b41 92000  //b51 90000
  float r1 = 500;

  int xmi1 = (int)(x1-r1);
  int xma1 = (int)(x1+r1);
  int ymi1 = (int)(y1-r1);
  int yma1 = (int)(y1+r1);

  int xmax1 = h2.GetXaxis()->FindBin(xma1);
  int xmin1 = h2.GetXaxis()->FindBin(xmi1);
  int ymax1 = h2.GetYaxis()->FindBin(yma1);
  int ymin1 = h2.GetYaxis()->FindBin(ymi1);

  float bkg1 = h2.Integral(xmin1,xmax1,ymin1,ymax1)/(xmax1-xmin1+1)/(ymax1-ymin1+1);
  *bkg = bkg1;
  // printf("Average background %f\n", *bkg);
}

void getMax(TH2F &h2, TObjArray &peaks, TObjArray &txt) {
  Int_t MaxBin = h2.GetMaximumBin();
  Int_t ix,iy,iz;
  int radius = 250;
  h2.GetBinXYZ(MaxBin, ix, iy, iz);
  float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
  float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
  TEllipse  *el = new TEllipse(x,y,radius,radius);
  el->SetFillStyle(0);
  peaks.Add(el);
  TText  *t = new TText(x,y+300,Form("%d",peaks.GetEntries()));
  //t->SetTextSize(0.8);
  txt.Add(t);
  int r0=4;
  for(int iix = ix-r0; iix<=ix+r0; iix++)
    for(int iiy = iy-r0; iiy<=iy+r0; iiy++)
      h2.SetBinContent(iix,iiy,0);
}

void get_peaks(TH2F &h2, TObjArray &peaks, TObjArray &txt, int npmax) {
  TH2F *h2new = (TH2F*)h2.Clone("get_peaks");
  for(int i=0; i<npmax; i++) getMax(*h2new, peaks, txt);
}

void drawEllipse(TObjArray &peaks, TObjArray &txt, int col) {
  int np = peaks.GetEntries();
  TText t(0,0,"a");
  for(int j=0; j<np; j++) {
    TEllipse *l = ((TEllipse*)(peaks.At(j)));
    l->SetLineColor(col);
    l->Draw();
    if(col!=kWhite) ((TText*)(txt.At(j)))->Draw();
  }
}

void getEntriesInEllipse(TH2F &h2, TEllipse &el, int *signif_bins, float *entries, float *bkg) {
  int nBinsX = h2.GetNbinsX();
  int nBinsY = h2.GetNbinsY();
  float x0 = el.GetX1();
  float y0 = el.GetY1();
  float r1 = el.GetR1();
  *signif_bins = 0;
  *entries = 0;

  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
      float x = h2.GetXaxis()->GetBinCenter(i);
      float y = h2.GetYaxis()->GetBinCenter(j);
      float distance = (pow(x-x0,2)+pow(y-y0,2))/pow(r1,2);
      if (distance <= 1) {
        int binContent = h2.GetBinContent(i, j);
        // float excess = (binContent - bkg)/TMath::Sqrt(bkg);
        // if (excess > 1) {
        *entries = *entries + binContent - *bkg;
          // *signif_bins = *signif_bins + 1;
        }
      }
    }
  // }
  // *entries = *entries - bkg;
  if (*entries<0) *entries = 0;
}

void count_bins(int npmax, TH2F &h2, TObjArray &peaks, int plate, TH1F *h_long[npmax]) {
  int np = peaks.GetEntries();
  int signif_bins;
  float bkg;
  float entries;
  float width;
  float bin_size = 50;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    eval_bkg(h2, &bkg);
    getEntriesInEllipse(h2, *el, &signif_bins, &entries, &bkg);
    // printf("Peak %i has %i significant bins and %i entries\n", i+1, signif_bins, (int)entries);
    // entries = entries - bkg;
    if (entries<0) entries = 0;
    h_long[i]->SetBinContent(plate, entries);
    width = TMath::Sqrt(signif_bins * bin_size * bin_size / TMath::Pi());
  }
}

void makePlots(TCanvas *c, int np, int npmax, TH1F* h_long, int *maxPeak, int *maxPlate) {
  int idx = (np%3) +1;
  if (idx==1) c->Clear("D");
  c->cd(idx)->SetGrid(1,0);
  *maxPeak = h_long->GetMaximum(),
  *maxPlate = h_long->GetMaximumBin();
  h_long->SetLineColor(1);
  h_long->SetLineWidth(2);
  h_long->Draw("hist");
  if(np == 2)                 c->Print("longitudinal_xz.pdf(","pdf");
  else if(np == npmax-1)      c->Print("longitudinal_xz.pdf)","pdf");
  else if(idx == 3 && np >2 ) c->Print("longitudinal_xz.pdf","pdf");
}

void findStart(TH1F* h_long, int *firstPlate, int *lastPlate) { //add second axis?
  h_long->Scale(1./h_long->Integral());
  *firstPlate = h_long->FindFirstBinAbove(0.01);
  *lastPlate = h_long->FindLastBinAbove(0.01);
}

void makeNtuple(int np, TH1F *h_long[np], TObjArray &peaks) {
  TFile *output = new TFile("peaks.root","RECREATE");
  TNtuple *ntuple = new TNtuple("showers","tagged showers","x:y:start:end:peak:maxplate");
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 1500);
  c2->Divide(1,3);
  
  float entries;
  float width;
  float bin_size = 50;
  int firstPlate, lastPlate, maxPeak, maxPlate;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    float x = el->GetX1();
    float y = el->GetY1();
    makePlots(c2, i, np, h_long[i], &maxPeak, &maxPlate);
    findStart(h_long[i], &firstPlate, &lastPlate);
    ntuple->Fill(x, y, firstPlate, lastPlate, maxPeak, maxPlate);
  }

  output->Write();
  output->Close();
}

void tag_basektrack() {
  int ntag = 50;
  TH2F *h2 = (TH2F *)(gDirectory->Get("XYseg"));
  TH3F *h3 = (TH3F *)(gDirectory->Get("XYPseg"));
  h2->Smooth();
  float bkg = set_limits(*h3,*h2);
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  c->SetGrid();
  gStyle->SetOptStat(0);
  h2->Draw("colz");
  c->Update();
  c->Print("sh.gif+180");
  
  TObjArray peaks;
  TObjArray txt;
  get_peaks(*h2,peaks,txt,ntag);
  drawEllipse(peaks,txt, kBlack);
  c->Update();
  c->Print("sh.gif+180");
  
  TH1F *h_long[ntag];
  for(int i=0; i<ntag; i++) {
    h_long[i] = new TH1F(Form("h_long_%i", i+1),Form("Cluster %i;plate;%%segments", i+1), 60, 1, 61);
  }

  const int nplates = 60; //change nplates for mc
  for(int p=1; p<=nplates; p++) { 
    h3->GetZaxis()->SetRange(p,p);
    TH2F *h = (TH2F*)(h3->Project3D("yx"));
    h->SetTitle(Form("Plate %d",p));
    h->Draw("colz");

    drawEllipse(peaks,txt, kBlack);
    count_bins(ntag, *h, peaks, p, &h_long[0]);
    c->Update();
    c->Print("sh.gif+12");
  }
  c->Print("sh.gif++");

  // makePlots(ntag, &h_long[0]);
  makeNtuple(ntag, &h_long[0], peaks);
}