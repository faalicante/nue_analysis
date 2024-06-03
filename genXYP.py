import ROOT
import fedrarootlogon
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cell", dest="cellID", required=False, default=None, type=int)
options = parser.parse_args()

brick = 21
nx = 18
ny = nx
xmin = 5000.
xmax = 185000.
ymin = 5000.
ymax = 185000.
emulsioncell = ROOT.EdbCell2()
emulsioncell.InitCell(nx,xmin,xmax,ny,ymin,ymax,1)
icell = options.cellID
ix=((icell % nx)) + 1
iy=((icell // ny)) + 1
xmincell = emulsioncell.X(ix-1)-emulsioncell.Xbin()
xmaxcell = emulsioncell.X(ix-1)+emulsioncell.Xbin()
ymincell = emulsioncell.Y(iy-1)-emulsioncell.Ybin()
ymaxcell = emulsioncell.Y(iy-1)+emulsioncell.Ybin()
bin_size = 50
xbin = int(3*(xmaxcell-xmincell)/bin_size)
ybin = int(3*(ymaxcell-ymincell)/bin_size)
overlap_fraction = 0.25

rootfile = ROOT.TFile(f"hist_XYP_{ix}_{iy}.root","RECREATE")
hXY = ROOT.TH2D(f"XYseg",f"XYseg;x[#mum];y[#mum]", xbin, xmincell, xmaxcell, ybin, ymincell, ymaxcell)
hXYP = ROOT.TH3D(f"XYPseg",f"XYPseg;x[#mum];y[#mum];p", xbin, xmincell, xmaxcell, ybin, ymincell, ymaxcell,60,1,60)
hTXTY = ROOT.TH2F("TXTYseg", "TXTYseg", 200, -0.1, 0.1, 200, -0.1, 0.1)

path = f"/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b{brick:06}"
#scanset
sspath = path                   ##CREATE FOLDERS
sproc = ROOT.EdbScanProc()
sproc.eProcDirClient=sspath
id = ROOT.EdbID(brick,0,0,0)
ss = sproc.ReadScanSet(id)
ss.Brick().SetID(brick)

npl = ss.eIDS.GetEntries()
cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
cutstring = cutstring + "&&" + f"TMath::Abs(s.eX-{emulsioncell.X(ix-1)}) < {emulsioncell.Xbin()*(1.+overlap_fraction)/2} && TMath::Abs(s.eY-{emulsioncell.Y(iy-1)}) < {emulsioncell.Ybin()*(1.+overlap_fraction)/2}"
cut = ROOT.TCut(cutstring)
cut.Print()
for i in range(npl):
  idplate = ss.GetID(i)
  nplate = idplate.ePlate
  plate = ss.GetPlate(idplate.ePlate)
  p = ROOT.EdbPattern()
  ect = ROOT.EdbCouplesTree()
  ect.InitCouplesTree("couples",path+f"/p{nplate:03}/{brick}.{nplate}.{ix}.{iy}.cp.root","READ")
  #addingcut
  ect.eCut = cut 
  cutlist = ect.InitCutList()
  if not ect.eTree: continue #no TTree in file
  if not cutlist:
    print("We have no entries, quitting!")
    continue
  nsegcut = cutlist.GetN()
  print("We have", nsegcut, "good couples in plate", nplate)
  
  # nseg = ect.eTree.GetEntries()
  
  for ientry in range(nsegcut):
    iseg = cutlist.GetEntry(ientry)
    ect.GetEntry(iseg)
    seg=ect.eS
    #setting z and affine transformation
    seg.SetZ(plate.Z())
    seg.SetPID(i)
    seg.Transform(plate.GetAffineXY())
    seg.Transform(plate.GetAffineTXYY())
    hXY.Fill(seg.X(), seg.Y())
    hXYP.Fill(seg.X(), seg.Y(),nplate)
    hTXTY.Fill(seg.TX(), seg.TY())
ect.Close()

rootfile.cd()
hXY.Write()
hXYP.Write()
hTXTY.Write()
rootfile.Close()



