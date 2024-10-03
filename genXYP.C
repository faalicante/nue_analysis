int genXYP(int cell){
  const int brick = 121;
  const int nx = 18;
  const int ny = nx;
  const float xmin = 5000.;
  const float xmax = 185000.;
  const float ymin = 5000.;
  const float ymax = 185000.;

  EdbCell2 *emulsioncell = new EdbCell2();
  emulsioncell->InitCell(nx,xmin,xmax,ny,ymin,ymax,1);
  const int ix=(cell % nx);
  const int iy=(cell / ny);
  const float xmincell = emulsioncell->X(ix)-emulsioncell->Xbin()/2;
  const float xmaxcell = emulsioncell->X(ix)+emulsioncell->Xbin()/2;
  const float ymincell = emulsioncell->Y(iy)-emulsioncell->Ybin()/2;
  const float ymaxcell = emulsioncell->Y(iy)+emulsioncell->Ybin()/2;
  const int bin_size = 50;
  const int xbin = (int)((xmaxcell-xmincell)/bin_size);
  const int ybin = (int)((ymaxcell-ymincell)/bin_size);
  const float overlap_fraction = 0.7;

  TFile *rootfile = new TFile(Form("hist_XYP_b%i_%i_%i.root", brick, ix, iy),"RECREATE");
  TH2D *hXY = new TH2D("XYseg","XYseg;x[#mum];y[#mum]", xbin, xmincell, xmaxcell, ybin, ymincell, ymaxcell);
  TNtuple *ntuple = new TNtuple("segments","Ntuple of segments","p:x:y:tx:ty:theta");
  TH2D *hXYs[57];
  for(int i=0; i<57; i++) {
    hXYs[57-i] = new TH2D(Form("XYseg_%i",57-i),Form("XYPseg_%i;x[#mum];y[#mum]",57-i), xbin, xmincell, xmaxcell, ybin, ymincell, ymaxcell);
  }
  TString path = Form("/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b%06i", brick);
  //scanset
  TString sspath = path+"/..";             
  EdbScanProc *sproc = new EdbScanProc();
  sproc->eProcDirClient=sspath;
  EdbID *id = new EdbID(brick,100,0,0);
  EdbScanSet *ss = sproc->ReadScanSet(*id);
  ss->Brick().SetID(brick);

  int npl = ss->eIDS.GetEntries();
  TString cutstring = "eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0";
  cutstring = cutstring + "&&" + Form("TMath::Abs(s.eX-%f) < %f && TMath::Abs(s.eY-%f) < %f", emulsioncell->X(ix), emulsioncell->Xbin()*(0.5+overlap_fraction), emulsioncell->Y(iy), emulsioncell->Ybin()*(0.5+overlap_fraction));
  TCut *cut = new TCut(cutstring);
  cut->Print();

  for (int i=0; i<npl; i++){
    EdbID *idplate = ss->GetID(i);
    int nplate = idplate->ePlate;
    EdbLayer *plate = ss->GetPlate(nplate);
    TString cp_file = path+Form("/p%03i/%i.%i.0.0.cp.root", nplate, brick, nplate);
    EdbCouplesTree *ect = new EdbCouplesTree();
    ect->InitCouplesTree("couples",cp_file,"READ");
    ect->eCut = *cut;
    if (!(ect->eTree)) return -1;
    TEventList *cutlist = ect->InitCutList();
    if (!cutlist){ 
      cout << "We have no entries, quitting!" << endl;
      return -1;
      }
    int nsegcut = cutlist->GetN();
    cout << "We have "<< nsegcut<< " good couples in plate "<< nplate << endl;
    
    int iseg = 0;
    for (int ientry=0; ientry < nsegcut; ientry++){
      iseg = cutlist->GetEntry(ientry);
      ect->GetEntry(iseg);
      EdbSegP *seg = ect->eS;
      seg->SetZ(plate->Z());
      seg->SetPID(i);
      seg->Transform(plate->GetAffineXY());
      EdbAffine2D *afftxty = plate->GetAffineTXTY();
      float tx = afftxty->A11()*seg->TX() + afftxty->A12()*seg->TY() + afftxty->B1();
      float ty = afftxty->A21()*seg->TX() + afftxty->A22()*seg->TY() + afftxty->B2();
      seg->SetTX(tx);
      seg->SetTY(ty);
      float sx = seg->X();
      float sy = seg->Y();
      if (TMath::Abs(sx-emulsioncell->X(ix)) < emulsioncell->Xbin()/2 && TMath::Abs(sy-emulsioncell->Y(iy)) < emulsioncell->Ybin()/2){
        hXY->Fill(sx, sy);
        hXYs[57-i]->Fill(sx, sy);
        ntuple->Fill(nplate, sx, sy, tx, ty, seg->Theta());
      }
    }
    rootfile->cd();
    hXYs[57-i]->Write();
    ect->Close();
  }

  rootfile->cd();
  hXY->Write();
  ntuple->Write();
  rootfile->Close();

  return 0;
}

