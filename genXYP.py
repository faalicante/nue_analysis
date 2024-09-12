import ROOT
import fedrarootlogon
import os.path

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", dest="fragment", required=False, default=None, type=int)
options = parser.parse_args()

brick = 121
nx = 18
ny = nx
xmin = 5000
xmax = 185000
ymin = 5000
ymax = 185000
# emulsioncell = ROOT.EdbCell2()
# emulsioncell.InitCell(nx,xmin,xmax,ny,ymin,ymax,1)
fragment = options.fragment
# ix=(icell % nx)
# iy=(icell // ny)
# xmincell = emulsioncell.X(ix)-emulsioncell.Xbin()
# xmaxcell = emulsioncell.X(ix)+emulsioncell.Xbin()
# ymincell = emulsioncell.Y(iy)-emulsioncell.Ybin()
# ymaxcell = emulsioncell.Y(iy)+emulsioncell.Ybin()
bin_size = 50
xbin = int((xmax-xmin)/bin_size)
ybin = int((ymax-ymin)/bin_size)
# overlap_fraction = 0.25

hXYs = {}
rootfile = ROOT.TFile(f"hist_XYP_run1b121_{fragment}.root","RECREATE")
hXY = ROOT.TH2D(f"XYseg",f"XYseg;x[#mum];y[#mum]", xbin, xmin, xmax, ybin, ymin, ymax)
hTXTY = ROOT.TH2F("TXTYseg", "TXTYseg", 200, -0.1, 0.1, 200, -0.1, 0.1)
for p in range(1, 58):
  hXYs[p] = ROOT.TH2D(f"XYseg_{p}",f"XYPseg_{p};x[#mum];y[#mum]", xbin, xmin, xmax, ybin, ymin, ymax)

path = f"/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b{brick:06}"
#scanset
sspath = path+f"/.."                  
sproc = ROOT.EdbScanProc()
sproc.eProcDirClient=sspath
id = ROOT.EdbID(brick,100,0,0)
ss = sproc.ReadScanSet(id)
ss.Brick().SetID(brick)

npl = ss.eIDS.GetEntries()
cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
# cutstring = cutstring + "&&" + f"TMath::Abs(s.eX-{emulsioncell.X(ix)}) < {emulsioncell.Xbin()*(1.+overlap_fraction)/2} && TMath::Abs(s.eY-{emulsioncell.Y(iy)}) < {emulsioncell.Ybin()*(1.+overlap_fraction)/2}"
cut = ROOT.TCut(cutstring)
cut.Print()
for i in range(npl):
  idplate = ss.GetID(i)
  nplate = idplate.ePlate
  plate = ss.GetPlate(idplate.ePlate)

  cp_file = path+f"/p{nplate:03}/{brick}.{nplate}.0.0.{fragment}.cp.root"
  if not os.path.isfile(cp_file): 
    print(f"Missin file /p{nplate:03}/{brick}.{nplate}.0.0.{fragment}.cp.root")
    continue

  p = ROOT.EdbPattern()
  ect = ROOT.EdbCouplesTree()
  ect.InitCouplesTree("couples",cp_file,"READ")
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
    afftxty = plate.GetAffineTXTY()
    tx = afftxty.A11()*seg.TX() + afftxty.A12()*seg.TY() + afftxty.B1()
    ty = afftxty.A21()*seg.TX() + afftxty.A22()*seg.TY() + afftxty.B2()
    seg.SetTX(tx)
    seg.SetTY(ty)
    hXY.Fill(seg.X(), seg.Y())
    hTXTY.Fill(seg.TX(), seg.TY())
    hXYs[nplate].Fill(seg.X(), seg.Y())
ect.Close()

rootfile.cd()
hXY.Write()
hTXTY.Write()
for p in range(1,58):
  hXYs[p].Write()
rootfile.Close()



