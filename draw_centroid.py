import ROOT
import fedrarootlogon
from math import pi
import numpy as np

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-b",dest="brick", type=int, required=True)
options = parser.parse_args()
brick = options.brick


def isShower(peak, seg):
  x0 = peak.x
  y0 = peak.y
  r0 = 300
  distance = ((x-x0)**2 + (y-y0)**2) / r0**2
  return distance <= 1


brick=121
path = f'/eos/user/f/falicant/Simulations_sndlhc/numucc_withcrisfiles_25_July_2022/b{brick:06}'
shower_file = ROOT.TFile.Open(path+f"/peaks.root")

peaks = shower_file.showers

histofile = ROOT.TFile(path+"/hist_couples_aligned.root","RECREATE")

cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
cut = ROOT.TCut(cutstring)
cut.Print()
for peak in peaks:
  tag = peak.tag
  print(f"Processing peak {tag}")
  startingPlate = peak.start
  endingPlate = peak.end
  for plate in range(1,61):
    if plate < startingPlate or plate > endingPlate: continue
    ect = ROOT.EdbCouplesTree()
    ect.InitCouplesTree("couples",path+f"/p{plate:03}/{brick}.{plate}.0.0.cp.root","READ")
    ect.eCut = cut
    cutlist = ect.InitCutList()
    nsegcut = cutlist.GetN()
    xmin=42400 - 2000
    xmax=42400 + 2000
    ymin=13000 - 2000
    ymax=13000 + 2000

    bin_size = 50
    xbin = int((xmax-xmin)/bin_size)
    ybin = int((ymax-ymin)/bin_size)
    hxy = ROOT.TH2D(f"h_{plate}",f"couples_plate_{plate};x[#mum];y[#mum]", xbin, xmin, xmax, ybin, ymin, ymax)
    for ientry in range(nsegcut):
      iseg = cutlist.GetEntry(ientry)
      ect.GetEntry(iseg)
      seg=ect.eS
      x = seg.X()
      y = seg.Y()
      tx = seg.TX()
      ty = seg.TY()
      if isInCluster(peak, x, y):
        hxy.Fill(x, y)
    ect.Close()
    hxy.Write()
histofile.Close()
