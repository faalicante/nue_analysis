import ROOT
import fedrarootlogon
import os.path

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cell", dest="cell", required=False, default=None, type=int)
options = parser.parse_args()

brick = 121
nx = 18
ny = nx
xmin = 5000.
xmax = 185000.
ymin = 5000.
ymax = 185000.
emulsioncell = ROOT.EdbCell2()
emulsioncell.InitCell(nx,xmin,xmax,ny,ymin,ymax,1)
cell = options.cell
ix=(cell % nx)
iy=(cell // ny)
xmincell = emulsioncell.X(ix)-emulsioncell.Xbin()/2
xmaxcell = emulsioncell.X(ix)+emulsioncell.Xbin()/2
ymincell = emulsioncell.Y(iy)-emulsioncell.Ybin()/2
ymaxcell = emulsioncell.Y(iy)+emulsioncell.Ybin()/2
bin_size = 50
xbin = int((xmaxcell-xmincell)/bin_size)
ybin = int((ymaxcell-ymincell)/bin_size)
overlap_fraction = 0.6

hXYs = {}
rootfile = ROOT.TFile(f"hist_XY_b{brick}_{cell}.root","RECREATE")
hXY = ROOT.TH2D(f"XYseg",f"XYseg;x[#mum];y[#mum]", xbin, xmincell, xmaxcell, ybin, ymincell, ymaxcell)
ntuple = ROOT.TNtuple("segments","Ntuple of segments","p:x:y:tx:ty:theta")

# hTXTY = ROOT.TH2F("TXTYseg", "TXTYseg", 200, -0.1, 0.1, 200, -0.1, 0.1)

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
cutstring = cutstring + "&&" + f"TMath::Abs(s.eX-{emulsioncell.X(ix)}) < {emulsioncell.Xbin()*(0.5+overlap_fraction)} && TMath::Abs(s.eY-{emulsioncell.Y(iy)}) < {emulsioncell.Ybin()*(0.5+overlap_fraction)}"
cut = ROOT.TCut(cutstring)
cut.Print()
for i in range(npl):
  idplate = ss.GetID(i)
  nplate = idplate.ePlate
  plate = ss.GetPlate(idplate.ePlate)

  cp_file = path+f"/p{nplate:03}/{brick}.{nplate}.0.0.cp.root"
  if not os.path.isfile(cp_file): 
    print(f"Missin file /p{nplate:03}/{brick}.{nplate}.0.0.cp.root")
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
  
  hXYs[nplate] = ROOT.TH2D(f"XYseg_{nplate}",f"XYseg_{nplate};x[#mum];y[#mum]", xbin, xmincell, xmaxcell, ybin, ymincell, ymaxcell)
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
    sx = seg.X()
    sy = seg.Y()
    if ROOT.TMath.Abs(sx-emulsioncell.X(ix)) < emulsioncell.Xbin()/2 and ROOT.TMath.Abs(sy-emulsioncell.Y(iy)) < emulsioncell.Ybin()/2:
      hXY.Fill(sx, sy)
      hXYs[nplate].Fill(sx, sy)
      ntuple.Fill(nplate, sx, sy, tx, ty, seg.Theta())
  rootfile.cd()
  hXYs[nplate].Write()
ect.Close()

rootfile.cd()
hXY.Write()
ntuple.Write()
rootfile.Close()



