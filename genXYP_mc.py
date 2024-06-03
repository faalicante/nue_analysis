import ROOT
import fedrarootlogon

brick = 21
xmin = 289000 - 5000
xmax = 299000 + 5000
ymin = 84000 - 5000
ymax = 94000 + 5000
bin_size = 50
xbin = int(3*(xmax-xmin)/bin_size)
ybin = int(3*(ymax-ymin)/bin_size)

rootfile = ROOT.TFile(f"hist_XYP_muon.root","RECREATE")
hXY = ROOT.TH2D(f"XYseg",f"XYseg;x[#mum];y[#mum]", xbin, xmin, xmax, ybin, ymin, ymax)
hXYP = ROOT.TH3D(f"XYPseg",f"XYPseg;x[#mum];y[#mum];p", xbin, xmin, xmax, ybin, ymin, ymax, 60, 1, 60)
hTXTY = ROOT.TH2F("TXTYseg", "TXTYseg", 200, -0.1, 0.1, 200, -0.1, 0.1)

cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
cutstring = cutstring + "&&" + f"TMath::Abs(s.eX - {(xmax+xmin)/2}) < {(xmax-xmin)/2} && TMath::Abs(s.eY - {(ymax+ymin)/2}) < {(ymax-ymin)/2}"
cut = ROOT.TCut(cutstring)
cut.Print()
for plate in range(1,61):
  ect = ROOT.EdbCouplesTree()
  ect.InitCouplesTree("couples",f"p{plate:03}/{brick}.{plate}.0.0.cp.root","READ")
  #addingcut
  ect.eCut = cut 
  cutlist = ect.InitCutList()
  if not ect.eTree: continue #no TTree in file
  if not cutlist:
    print("We have no entries, quitting!")
    continue
  nsegcut = cutlist.GetN()
  print("We have", nsegcut, "good couples in plate", plate)
  for ientry in range(nsegcut):
    iseg = cutlist.GetEntry(ientry)
    ect.GetEntry(iseg)
    seg=ect.eS
    hXY.Fill(seg.X(), seg.Y())
    hXYP.Fill(seg.X(), seg.Y(), plate)
    hTXTY.Fill(seg.TX(), seg.TY())
ect.Close()

rootfile.cd()
hXY.Write()
hXYP.Write()
hTXTY.Write()
rootfile.Close()



