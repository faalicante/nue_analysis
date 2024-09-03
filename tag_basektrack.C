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

void eval_bkg(TH2F &h2, TEllipse &el, float *bkg) { //make it two control regions next to the circle
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

void getMax(TH2F &h2, TObjArray &peaks, TObjArray &txt) {
  Int_t MaxBin = h2.GetMaximumBin();
  Int_t ix,iy,iz;
  int radius = 200;
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
  float bin_size = 50./3;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    eval_bkg(h2, *el, &bkg);
    getEntriesInEllipse(h2, *el, &signif_bins, &entries, &bkg);
    // printf("Peak %i has %i significant bins and %i entries\n", i+1, signif_bins, (int)entries);
    // entries = entries - bkg;
    if (entries<0) entries = 0;
    h_long[i]->SetBinContent(plate, entries);
    width = TMath::Sqrt(signif_bins * bin_size * bin_size / TMath::Pi());
  }
}

void makePlots(int npmax, TH1F* h_long[npmax]) {
  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 600);
  TLegend *leg = new TLegend(0.75, 0.7, 0.87, 0.87);
  leg->SetTextFont(53);
  leg->SetBorderSize(0);
  float max_long = 0;

  for(int i=0; i<npmax; i++) {
    leg->AddEntry(h_long[i], Form("Peak %i",i), "L2");
    h_long[i]->SetLineColor(i+1);
    h_long[i]->SetLineWidth(2);
    h_long[i]->Scale(1./h_long[i]->Integral());
    if (h_long[i]->GetMaximum() > max_long) {
      max_long = h_long[i]->GetMaximum();
      h_long[0]->GetYaxis()->SetRangeUser(0, max_long*1.2);
    }
    h_long[i]->Draw("same&&hist");
    c2->Update();

  }
  c2->Print("longitudinal_xz.png");
}

void findStart(TH1F* h_long, int *firstPlate, int *lastPlate) {
 *firstPlate = h_long->FindFirstBinAbove(0.005);
 *lastPlate = h_long->FindLastBinAbove(0.01);
}

void makeNtuple(int np, TH1F *h_long[np], TObjArray &peaks) {
  TFile *output = new TFile("peaks.root","RECREATE");
  TNtuple *ntuple = new TNtuple("showers","tagged showers","x:y:start:end");
  int signif_bins;
  
  float entries;
  float width;
  float bin_size = 50;
  int firstPlate, lastPlate;

  for(int i=0; i<np; i++) {
    TEllipse *el = ((TEllipse*)(peaks.At(i)));
    float x = el->GetX1();
    float y = el->GetY1();
    findStart(h_long[i], &firstPlate, &lastPlate);
    ntuple->Fill(x, y, firstPlate, lastPlate);
  }

  output->Write();
  output->Close();
}

void tag_basektrack() {
  int ntag = 5;
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
    h_long[i] = new TH1F(Form("h_long_%i", i+1),"Longitudinal entries;plate", 60, 1, 61);
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

  makePlots(ntag, &h_long[0]);
  makeNtuple(ntag, &h_long[0], peaks);
}