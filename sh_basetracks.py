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

vecs = {}
def append_list_to_dict(key, values_list):
    if key in vecs:
        vecs[key].append(values_list)
    else:
        vecs[key] = [values_list]

xmin_muon = 289000 - 2000
xmax_muon = 299000 + 2000
ymin_muon = 84000 - 2000
ymax_muon = 94000 + 2000
bin_size = 50
xbin = int((xmax_muon-xmin_muon)/bin_size)
ybin = int((ymax_muon-ymin_muon)/bin_size)


prepath = f'/eos/user/f/falicant/Simulations_sndlhc/muon1E5_simsndlhc'
rootfile = ROOT.TFile(prepath+f"/ssh_theta_{brick}.root","RECREATE")
ntuple = ROOT.TNtuple("couples", "Tree of couples","p:x:y:tx:ty:theta:scatt")

shower_file = ROOT.TFile.Open(prepath+"/muon_showers.root")
cbmsim = shower_file.cbmsim
shower_list = []
h = {}
for event in cbmsim:
  if event.count < 2000: continue
  shower_list.append(event.event)  
  h[f'xy_{event.event}'] = ROOT.TH2D(f'xy_{event.event}', f'xy_{event.event}; x [#mum]; y [#mum]', xbin, xmin_muon, xmax_muon, ybin, ymin_muon, ymax_muon)
  h[f'theta_{event.event}'] = ROOT.TH2D(f'theta_{event.event}', f'theta_{event.event}; plate; theta [#mrad]', 60, 1, 61, 3140, 0, pi*100)
h['xy'] = ROOT.TH2D('xy', 'xy; x [#mum]; y [#mum]', xbin, xmin_muon, xmax_muon, ybin, ymin_muon, ymax_muon)
h['theta_mean'] = ROOT.TH1D('theta_mean', 'theta_mean; theta [#mrad]', 3140, 0, pi*100)
h['theta_std'] = ROOT.TH1D('theta_std', 'theta_std; theta [#mrad]', 3140, 0, pi*100)
h['spherocity'] = ROOT.TH1D('spherocity', 'spherocity', 100, 0, 1)
h['direction'] = ROOT.TH1D('direction', 'direction; phi', 2*3140, 0, 2*pi)
print('Selected muon shower events:', len(shower_list))

cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
cut = ROOT.TCut(cutstring)
# cut_muon.Print()

brick_path = prepath +f'/b{brick:06}'
for plate in range(1,61):
  ect = ROOT.EdbCouplesTree()
  ect.InitCouplesTree("couples",brick_path+f"/p{plate:03}/{brick}.{plate}.0.0.cp.root","READ")
  #addingcut
  ect.eCut = cut
  cutlist = ect.InitCutList()
  if not ect.eTree: continue #no TTree in file
  if not cutlist:
    print("We have no entries, quitting!")
    continue
  nsegcut = cutlist.GetN()
  print(f"Processing brick {brick} plate {plate}")
  for ientry in range(nsegcut):
    iseg = cutlist.GetEntry(ientry)
    ect.GetEntry(iseg)
    seg=ect.eS
    event = seg.MCEvt()
    if event not in shower_list: continue
    sx = seg.X()
    sy = seg.Y()
    stx = seg.TX()
    sty = seg.TY()
    append_list_to_dict(event, [stx, sty])
    stheta =seg.Theta()*100
    ntuple.Fill(plate, sx, sy, stx, sty , stheta)
    h[f'xy_{event}'].Fill(sx,sy)
    h[f'theta_{event}'].Fill(plate, stheta)

print(len(vecs), "event in this brick out of", len(shower_list)) 
ect.Close()
rootfile.cd()
for event in vecs:
  h[f'xy_{event}'].Write()
  h[f'theta_{event}'].Write()
  mean = h[f'theta_{event}'].GetMean(2)
  mean_err = h[f'theta_{event}'].GetMeanError(2)
  std = h[f'theta_{event}'].GetStdDev(2)
  std_err = h[f'theta_{event}'].GetStdDevError(2)
  h['theta_mean'].Fill(mean)
  h['theta_std'].Fill(std)
  sph, direction = spherocity(vecs[event])
  h['spherocity'].Fill(sph)
  h['direction'].Fill(direction[0])
h['theta_mean'].Write()
h['theta_std'].Write()
ntuple.Write()
rootfile.Close()

##issue: same event in multiple bricks
## do it by cluster tag