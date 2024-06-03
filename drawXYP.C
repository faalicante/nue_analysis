float set_range(TH3F &h3, TH2F &h2, float thr){
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

float set_limits(TH3F &h3, TH2F &h2){
  float meanbin = h2.Integral()/h2.GetNbinsX()/h2.GetNbinsY();
  float newmean = set_range(h3, h2, meanbin/2);
  float lastmean = set_range(h3, h2, newmean/2);
  printf("meanbin = %f    newmean = %f   lastmean = %f\n",meanbin, newmean, lastmean);
  return(lastmean)
}

void getMax(TH2F &h2, TObjArray &a, TObjArray &txt){
   Int_t MaxBin = h2.GetMaximumBin();
   Int_t ix,iy,iz;
   h2.GetBinXYZ(MaxBin, ix, iy, iz);
   float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
   float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
   TEllipse  *el = new TEllipse(x,y,300,300);
   el->SetFillStyle(0);
   a.Add(el);
   TText  *t = new TText(x,y+300,Form("%d",a.GetEntries()));
   //t->SetTextSize(0.8);
   txt.Add(t);
   //printf("%f %f \n",x,y);
   int r0=4;
   for(int iix = ix-r0; iix<=ix+r0; iix++)
     for(int iiy = iy-r0; iiy<=iy+r0; iiy++)
       h2.SetBinContent(iix,iiy,0);
}

void get_peaks(TH2F &h2, TObjArray &peaks, TObjArray &txt, int npmax){
  TH2F *h2new = (TH2F*)h2.Clone("get_peaks");
  for(int i=0; i<npmax; i++) getMax(*h2new, peaks, txt);
}

void drawEllipse(TObjArray &peaks, TObjArray &txt, int col){
  int np = peaks.GetEntries();
  TText t(0,0,"a");
  for(int j=0; j<np; j++) {
    TEllipse *l = ((TEllipse*)(peaks.At(j)));
    l->SetLineColor(col);
    l->Draw();
    if(col!=kWhite) ((TText*)(txt.At(j)))->Draw();
  }
}

void drawXYP(){
  TH2F *h2 = (TH2F *)(gDirectory->Get("XYseg"));
  TH3F *h3 = (TH3F *)(gDirectory->Get("XYPseg"));
  h2->Smooth();
  float bkg = set_limits(*h3,*h2);
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  c->SetGrid();
  gStyle->SetOptStat("n");
  h2->Draw("colz");
  c->Update();
  c->Print("sh.gif+180");
  
  TObjArray peaks;
  TObjArray txt;
  get_peaks(*h2,peaks,txt,5);
  drawEllipse(peaks,txt, kBlack);
  c->Update();
  c->Print("sh.gif+180");

  for(int i=1; i<=60; i++) { //change nplates for mc
    h3->GetZaxis()->SetRange(i,i);
    TH2F *h = (TH2F*)(h3->Project3D("yx"));
    h->SetTitle(Form("plate %d",i));
    h->Smooth();
    h->Draw("colz");
    drawEllipse(peaks,txt, kBlack);

    TObjArray locpeaks;
    TObjArray loctxt;
    get_peaks(*h,locpeaks,loctxt,10);
    drawEllipse(locpeaks,loctxt, kWhite);

    c->Update();
    c->Print("sh.gif+12");
  }
  c->Print("sh.gif++");
}

