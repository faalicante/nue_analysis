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


float eval_bkg(TH2F &h2) {
  int xmi = 291000;
  int xma = 297000;
  int ymi = 86000;
  int yma = 93000;
  // int xmi = 44000;
  // int xma = 57000;
  // int ymi = 42000;
  // int yma = 56000;
  // h2.GetXaxis()->SetRangeUser(xmi,xma);
  // h2.GetYaxis()->SetRangeUser(ymi,yma);
  int xmax = h2.GetXaxis()->FindBin(xma);
  int xmin = h2.GetXaxis()->FindBin(xmi);
  int ymax = h2.GetYaxis()->FindBin(yma);
  int ymin = h2.GetYaxis()->FindBin(ymi);
  float bkg = h2.Integral(xmin,xmax,ymin,ymax)/(xmax-xmin+1)/(ymax-ymin+1);
  printf("Average background %f\n", bkg);
  return bkg;
}

/*
float eval_bkg(TH2F &h2) {
  int xmi = 290000;
  int xma = 298000;
  int ymi = 85000;
  int yma = 94000;
  h2.GetXaxis()->SetRangeUser(xmi,xma);
  h2.GetYaxis()->SetRangeUser(ymi,yma);
  int binmax = h2.GetMaximum();
  int nBinsX = h2.GetNbinsX();
  int nBinsY = h2.GetNbinsY();
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  TH1F* hmode = new TH1F("hmode", "hmode", binmax, 1, binmax+1);
  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
        hmode->Fill(h2.GetBinContent(i,j));
    }
  }
  int modebin = hmode->GetMaximumBin();
  float bkg = hmode->GetXaxis()->GetBinLowEdge(modebin);
  printf("Mode background %f\n", bkg);
  hmode->Draw();
  c1->Draw();
  return bkg;
}
*/

void getMax(TH2F &h2, TObjArray &peaks, TObjArray &txt) {
  Int_t MaxBin = h2.GetMaximumBin();
  Int_t ix,iy,iz;
  h2.GetBinXYZ(MaxBin, ix, iy, iz);
  float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
  float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
  TEllipse  *el = new TEllipse(x,y,300,300);
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

void getEntriesInEllipse(TH2F &h2, TEllipse &el, float bkg, int *signif_bins, float *entries) {
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
        float excess = (binContent - bkg)/TMath::Sqrt(bkg);
        if (excess > 1) {
          *entries = *entries + binContent;
          *signif_bins = *signif_bins + 1;
        }
      }
    }
  }
  *entries = *entries - bkg;
  if (*entries<0) *entries = 0;
}

void count_bins(TH2F &h2, TObjArray &peaks, float bkg) {
  int np = peaks.GetEntries();
  int signif_bins;
  float entries;

  for(int j=0; j<np; j++) {
    TEllipse *el = ((TEllipse*)(peaks.At(j)));
    getEntriesInEllipse(h2, *el, bkg, &signif_bins, &entries);
    printf("Peak %i has %i significant bins and %i entries\n", j+1, signif_bins, (int)entries);
  }
}

void count_bins(TH2F &h2, TObjArray &peaks, float bkg, int plate, TH1F *h_prof[5], TH1F *h_signif[5], TH1F *h_width[5]) {
  int np = peaks.GetEntries();
  int signif_bins;
  float entries;
  float width;
  float bin_size = 50./3;

  for(int j=0; j<np; j++) {
    TEllipse *el = ((TEllipse*)(peaks.At(j)));
    getEntriesInEllipse(h2, *el, bkg, &signif_bins, &entries);
    printf("Peak %i has %i significant bins and %i entries\n", j+1, signif_bins, (int)entries);
    h_prof[j]->SetBinContent(plate, entries);
    h_signif[j]->SetBinContent(plate, signif_bins);
    width = TMath::Sqrt(signif_bins * bin_size * bin_size / TMath::Pi());
    h_width[j]->SetBinContent(plate, width);
  }
}

void makePlots(int npmax, TH1F* h_prof[5], TH1F* h_signif[5], TH1F* h_width[5]) {
  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 600);
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 600);
  TCanvas *c4 = new TCanvas("c4", "c4", 1000, 600);
  TLegend *leg = new TLegend(0.75, 0.7, 0.87, 0.87);
  leg->SetTextFont(53);
  leg->SetBorderSize(0);
  int line_color[5] = {1, 2, 3, 4, 6};
  float max_prof = 0;
  float max_signif = 0;
  float max_width = 0;

  for(int j=0; j<npmax; j++) {
    leg->AddEntry(h_prof[j], Form("Peak %i",j), "L2");
    c2->cd();
    if (h_prof[j]->GetMaximum() > max_prof) {
      max_prof = h_prof[j]->GetMaximum();
      h_prof[0]->GetYaxis()->SetRangeUser(0, max_prof*1.2);
    }
    h_prof[j]->SetLineColor(line_color[j]);
    h_prof[j]->SetLineWidth(2);
    h_prof[j]->Draw("same");
    c3->cd();
    if (h_signif[j]->GetMaximum() > max_signif) {
      max_signif = h_signif[j]->GetMaximum();
      h_signif[0]->GetYaxis()->SetRangeUser(0, max_signif*1.2);
    }
    h_signif[j]->SetLineColor(line_color[j]);
    h_signif[j]->SetLineWidth(2);
    h_signif[j]->Draw("same");
    c4->cd();
    if (h_width[j]->GetMaximum() > max_width) {
      max_width = h_width[j]->GetMaximum();
      h_width[0]->GetYaxis()->SetRangeUser(0, max_width*1.2);
    }
    h_width[j]->SetLineColor(line_color[j]);
    h_width[j]->SetLineWidth(2);
    h_width[j]->Draw("same");
    c2->Update();
    c3->Update();
    c4->Update();
  }
  c2->Print("prof.png");
  c3->Print("signif.png");
  c4->Print("width.png");
}

void drawXYP() {
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
  get_peaks(*h2,peaks,txt,5);
  drawEllipse(peaks,txt, kBlack);
  c->Update();
  c->Print("sh.gif+180");
  // float bkg = h2->Integral()/h2->GetNbinsX()/h2->GetNbinsY();
  // float bkg = eval_bkg(*h2);
  count_bins(*h2, peaks, bkg);
  
  int npmax = peaks.GetEntries();
  TH1F *h_prof[npmax], *h_signif[npmax], *h_width[npmax];
  for(int j=0; j<npmax; j++) {
    h_prof[j] = new TH1F(Form("h_prof_%i", j+1),"Profile entries;plate", 60, 1, 61);
    h_signif[j] = new TH1F(Form("h_signif_%i", j+1),"Profile significant bins;plate", 60, 1, 61);
    h_width[j] = new TH1F(Form("h_width_%i", j+1),"Profile width;plate", 60, 1, 61);
  }

  const int nplates = 60; //change nplates for mc
  for(int i=1; i<=nplates; i++) { 
    h3->GetZaxis()->SetRange(i,i);
    TH2F *h = (TH2F*)(h3->Project3D("yx"));
    h->SetTitle(Form("Plate %d",i));
    h->Smooth();
    //h->Rebin2D(3,3);
    bkg = eval_bkg(*h);
    // bkg = h->Integral()/h->GetNbinsX()/h->GetNbinsY();
    // printf("Average background %f\n", bkg);

    h->Draw("colz");

    drawEllipse(peaks,txt, kBlack);

    TObjArray locpeaks;
    TObjArray loctxt;
    get_peaks(*h,locpeaks,loctxt,10);
    drawEllipse(locpeaks,loctxt, kWhite);

    count_bins(*h, peaks, bkg, i, &h_prof[0], &h_signif[0], &h_width[0]);
    c->Update();
    c->Print("sh.gif+12");
  }
  c->Print("sh.gif++");
  makePlots(npmax, &h_prof[0], &h_signif[0], &h_width[0]);

  // TFile *output_file = new TFile("prof.root","RECREATE");
  // for(int j=0; j<npmax; j++) {
  //   h_prof[j]->Write();
  //   h_signif[j]->Write();
  //   h_width[j]->Write();
  // }
  // output_file->Write();
  // output_file->Close();
}

