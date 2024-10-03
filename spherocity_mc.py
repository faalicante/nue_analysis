import ROOT
import fedrarootlogon
from math import pi
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-b",dest="brick", type=int, required=True)
options = parser.parse_args()
brick = options.brick

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

def isInCluster(peak, x, y):
  x0 = peak.x
  y0 = peak.y
  r0 = 300
  distance = ((x-x0)**2 + (y-y0)**2) / r0**2
  return distance <= 1

# path = f'/eos/user/f/falicant/Simulations_sndlhc/muon1E5_simsndlhc/b{brick:06}'
# path = f'/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b{brick:06}'
path = f'/eos/user/f/falicant/Simulations_sndlhc/numucc_withcrisfiles_25_July_2022/b{brick:06}'

shower_file = ROOT.TFile.Open(path+f"/peaks.root")
peaks = shower_file.showers
h = {}
rootfile = ROOT.TFile(path+f"/sh_angle_{brick}.root","RECREATE")
ntuple = ROOT.TNtuple("couples", "Tree of couples","brick:peak:max:maxpeak:start:end:nseg:theta:spherocity:dir:event:rank")
h['theta_mean'] = ROOT.TH1D('theta_mean', 'theta_mean; theta [#mrad]', 150, 0, 150)
h['spherocity'] = ROOT.TH1D('spherocity', 'spherocity', 100, 0, 1)
h['direction'] = ROOT.TH1D('direction', 'direction; phi', 150, 0, 150)

cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
cut = ROOT.TCut(cutstring)
cut.Print()
for ipeak, peak in enumerate(peaks):
  DictMCEvt = {}
  print(f"Processing peak {ipeak+1}")
  vecs = []
  startingPlate = peak.start
  if startingPlate < 5: continue
  endingPlate = peak.end
  h[f'theta_{ipeak+1}'] = ROOT.TH1D(f'theta_{ipeak+1}', f'theta_{ipeak+1}; theta [#mrad]',150, 0, 150)
  for plate in range(1,61):
    if plate < startingPlate or plate > endingPlate: continue
    ect = ROOT.EdbCouplesTree()
    ect.InitCouplesTree("couples",path+f"/p{plate:03}/{brick}.{plate}.0.0.cp.root","READ")
    ect.eCut = cut
    cutlist = ect.InitCutList()
    nsegcut = cutlist.GetN()
    for ientry in range(nsegcut):
      iseg = cutlist.GetEntry(ientry)
      ect.GetEntry(iseg)
      seg=ect.eS
      x = seg.X()
      y = seg.Y()
      if isInCluster(peak, x, y):
        tx = seg.TX()
        ty = seg.TY()
        vecs.append([tx, ty])
        theta =seg.Theta()*1000
        h[f'theta_{ipeak+1}'].Fill(theta)
        sEvt = seg.MCEvt()
        if sEvt in DictMCEvt:
            DictMCEvt[sEvt] += 1
        else:
            DictMCEvt[sEvt] = 1
    ect.Close()
  eventID = max(DictMCEvt, key=DictMCEvt.get)
  sph, direction = spherocity(vecs)
  mean = h[f'theta_{ipeak+1}'].GetMean()
  h['theta_mean'].Fill(mean)
  h['spherocity'].Fill(sph)
  h['direction'].Fill(direction[0])
  ntuple.Fill(brick, ipeak+1, peak.peak, peak.maxplate, startingPlate, endingPlate, peak.nseg, mean, sph, direction[0],eventID, peak.rankbin)

rootfile.cd()
h['theta_mean'].Write()
h['spherocity'].Write()
h['direction'].Write()
ntuple.Write()
rootfile.Close()
