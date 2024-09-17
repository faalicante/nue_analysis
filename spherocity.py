##IMPORTANT###
##NEED TO SAVE THE NTUPLE PER PEAK TO GET SPHEROCITY EVENT
import ROOT
import fedrarootlogon
from math import pi
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-b",dest="brick", type=int, required=True)
parser.add_argument("-f",dest="fragment", type=int, required=True)
options = parser.parse_args()
brick = options.brick
fragment = options.fragment

# Function to minimise
def cross(angle, vecs):
    c = np.cross([np.cos(angle[0]), np.sin(angle[0])], vecs)
    ret = np.sum(np.abs(c))
    return ret

# Spherocity
def spherocity(vecs, initial_angle = 0.):
    opt = scipy.optimize.minimize(cross, initial_angle, vecs, method = "Nelder-Mead")
    num = opt.fun
    denom = np.sum(np.linalg.norm(vecs, axis = 1))
    direction = opt.x
    return np.pi**2/4.*(num/denom)**2, direction

vecs = {}
def append_list_to_vecs(key, values_list):
    if key in vecs:
        vecs[key].append(values_list)
    else:
        vecs[key] = [values_list]

def isInCluster(peak, seg):
  x = seg.X()
  y = seg.Y()
  x0 = peak.x
  y0 = peak.y
  r0 = 250
  distance = ((x-x0)**2 + (y-y0)**2) / r0**2
  return distance <= 1


# path = f'/eos/user/f/falicant/Simulations_sndlhc/muon1E5_simsndlhc/b{brick:06}'
# path = f'/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b{brick:06}'
path = '/eos/user/f/falicant/nue_search/R1B121/30'
path = '/afs/cern.ch/work/f/falicant/public/nue_analysis/30'
datapath = '/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b000121'

shower_file = ROOT.TFile.Open(path+f"/peaks_{fragment}.root")
peaks = shower_file.showers
h = {}

#scanset
sspath = datapath+f"/.."
sproc = ROOT.EdbScanProc()
sproc.eProcDirClient=sspath
id = ROOT.EdbID(brick,100,0,0)
ss = sproc.ReadScanSet(id)
ss.Brick().SetID(brick)
npl = ss.eIDS.GetEntries()

cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
cut = ROOT.TCut(cutstring)
# cut.Print()

rootfile = ROOT.TFile(path+f"/sh_angle_{fragment}.root","RECREATE")
ntuple = ROOT.TNtuple("couples", "Tree of couples","fragment:peak:max:maxpeak:start:end:theta:spherocity")
h['theta_mean'] = ROOT.TH1D('theta_mean', 'theta_mean; theta [#mrad]', 150, 0, 150)
h['spherocity'] = ROOT.TH1D('spherocity', 'spherocity', 100, 0, 1)
h['direction'] = ROOT.TH1D('direction', 'direction; phi', 150, 0, 150)

for ipeak, peak in enumerate(peaks):
  vecs = []
  #if ipeak>4: break
  # xmin = peak.x - 500
  # xmax = peak.x + 500
  # ymin = peak.y - 500
  # ymax = peak.y + 500
  # bin_size = 50
  # xbin = int((xmax-xmin)/bin_size)
  # ybin = int((ymax-ymin)/bin_size)
  # h[f'xy_{ipeak}'] = ROOT.TH2D(f'xy_{ipeak}', f'xy_{ipeak}; x [#mum]; y [#mum]', xbin, xmin, xmax, ybin, ymin, ymax)
  h[f'theta_{ipeak}'] = ROOT.TH1D(f'theta_{ipeak}', f'theta_{ipeak}; theta [#mrad]',150, 0, 150)

  startingPlate = peak.start
  if startingPlate < 5: continue
  endingPlate = peak.end
  for i in range(npl):
    idplate = ss.GetID(i)
    plate = idplate.ePlate
    nplate = ss.GetPlate(plate)
  # for plate in range(1,61):
    if plate < startingPlate or plate > endingPlate: continue
    ect = ROOT.EdbCouplesTree()
    ect.InitCouplesTree("couples",datapath+f"/p{plate:03}/{brick}.{plate}.0.0.{fragment}.cp.root","READ")
    ect.eCut = cut
    cutlist = ect.InitCutList()
    if not ect.eTree: continue
    if not cutlist:
      print("We have no entries, quitting!")
      continue
    nsegcut = cutlist.GetN()
    print(f"Processing brick {brick} peak {ipeak} plate {plate}")
    for ientry in range(nsegcut):
      iseg = cutlist.GetEntry(ientry)
      ect.GetEntry(iseg)
      seg=ect.eS
#      event = seg.MCEvt()
      seg.SetZ(nplate.Z())
      seg.SetPID(i)
      seg.Transform(nplate.GetAffineXY())

      if isInCluster(peak, seg):
        # sx = seg.X()
        # sy = seg.Y()
        stx = seg.TX()
        sty = seg.TY()
        # append_list_to_vecs(ipeak, [stx, sty])
        vecs.append([stx, sty])
        stheta =seg.Theta()*1000
        # ntuple.Fill(brick, ipeak, plate, sx, sy, stx, sty , stheta)
        # h[f'xy_{ipeak}'].Fill(sx,sy)
        h[f'theta_{ipeak}'].Fill(stheta)
    ect.Close()
  sph, direction = spherocity(vecs)
  mean = h[f'theta_{ipeak}'].GetMean(2)
  h['theta_mean'].Fill(mean)
  h['spherocity'].Fill(sph)
  h['direction'].Fill(direction[0])
  ntuple.Fill(fragment, ipeak, peak.peak, peak.maxplate, startingPlate, endingPlate, mean, sph)

rootfile.cd()
# for event in vecs:
  # h[f'xy_{event}'].Write()
  # h[f'theta_{event}'].Write()
  # sph, direction = spherocity(vecs[event])
  # ntuple.Fill(fragment,event,mean,sph)
h['theta_mean'].Write()
h['spherocity'].Write()
h['direction'].Write()
ntuple.Write()
rootfile.Close()
