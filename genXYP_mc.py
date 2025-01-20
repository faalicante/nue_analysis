import ROOT
import fedrarootlogon

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-b",dest="brick", type=int, required=True)
# parser.add_argument("--cell",dest="cell", type=int, required=True)
options = parser.parse_args()
brick = options.brick
# cell = options.cell
# cellx = cell%10
# celly = cell//10
# cellsize = 20000

xmin = 289000 - 1000
xmax = 299000 + 1000
ymin = 84000 - 1000
ymax = 94000 + 1000
# bricks1 = [11,21,31,41,51]
# bricks2 = [12,22,32,42,52]
# bricks3 = [13,23,33,43,53]
# bricks4 = [14,24,34,44,54]
# if brick in bricks1:
#   offsetx = 195000
#   offsety = 0
# elif brick in bricks2:
#   offsetx = 15
#   offsety = 0
# elif brick in bricks3:
#   offsetx = 195000
#   offsety = 195000
# elif brick in bricks4:
#   offsetx = 15
#   offsety = 195000

# xmin = offsetx
# xmax = xmin + 195000
# ymin = offsety
# ymax = ymin + 195000
bin_size = 50
xbin = int((xmax-xmin)/bin_size)
ybin = int((ymax-ymin)/bin_size)

hXYs = {}
# rootfile = ROOT.TFile(f"hist_XY_numu.root","RECREATE")
# rootfile = ROOT.TFile(f"hist_XY_nue.root","RECREATE")
rootfile = ROOT.TFile(f"hist_XY_muon.root","RECREATE")
hXY = ROOT.TH2D(f"XYseg",f"XYseg;x[#mum];y[#mum]", xbin, xmin, xmax, ybin, ymin, ymax)
hTXTY = ROOT.TH2F("TXTYseg", "TXTYseg", 200, -0.1, 0.1, 200, -0.1, 0.1)
for p in range(1, 61):
  hXYs[p] = ROOT.TH2D(f"XYseg_{p}",f"XYseg_{p};x[#mum];y[#mum]", xbin, xmin, xmax, ybin, ymin, ymax)

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
    hTXTY.Fill(seg.TX(), seg.TY())
    hXYs[plate].Fill(seg.X(), seg.Y())
ect.Close()

rootfile.cd()
hXY.Write()
hTXTY.Write()
for p in range(1,61):
  hXYs[p].Write()
rootfile.Close()



