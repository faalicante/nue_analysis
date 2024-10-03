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
# parser.add_argument("-b",dest="brick", type=int, required=True)
parser.add_argument("--cell",dest="cell", type=int, required=True)
options = parser.parse_args()
# brick = options.brick
cell = options.cell
ix=(cell % 18)
iy=(cell // 18)

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

# vecs = {}
# def append_list_to_vecs(key, values_list):
#     if key in vecs:
#         vecs[key].append(values_list)
#     else:
#         vecs[key] = [values_list]

def isInCluster(peak, x, y):
  x0 = peak.x
  y0 = peak.y
  r0 = 300
  distance = ((x-x0)**2 + (y-y0)**2) / r0**2
  return distance <= 1


# path = f'/eos/user/f/falicant/Simulations_sndlhc/muon1E5_simsndlhc/b{brick:06}'
# path = f'/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b{brick:06}'
path = '/eos/user/f/falicant/nue_search/R1B121/gen1'
# path = '/afs/cern.ch/work/f/falicant/public/nue_analysis/30'
# datapath = '/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b000121'

shower_file = ROOT.TFile.Open(path+f"/peaks2_{cell}.root")
peaks = shower_file.showers
h = {}

#scanset
# sspath = datapath+f"/.."
# sproc = ROOT.EdbScanProc()
# sproc.eProcDirClient=sspath
# id = ROOT.EdbID(brick,100,0,0)
# ss = sproc.ReadScanSet(id)
# ss.Brick().SetID(brick)
# npl = ss.eIDS.GetEntries()

# cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
# cut = ROOT.TCut(cutstring)
# cut.Print()

cp_file = ROOT.TFile.Open(path+f"/hist/hist_XYP_b121_{cell}.root")
segments = cp_file.segments
rootfile = ROOT.TFile(path+f"/sh_angle_{cell}.root","RECREATE")
ntuple = ROOT.TNtuple("couples", "Tree of couples","cell:cellx:celly:peak:max:maxpeak:start:end:nseg:theta:spherocity:dir:rank")
h['theta_mean'] = ROOT.TH1D('theta_mean', 'theta_mean; theta [#mrad]', 150, 0, 150)
h['spherocity'] = ROOT.TH1D('spherocity', 'spherocity', 100, 0, 1)
h['direction'] = ROOT.TH1D('direction', 'direction; phi', 150, 0, 150)

for ipeak, peak in enumerate(peaks):
  print(f"Processing peak {ipeak+1}")
  vecs = []
  #if ipeak+1>4: break
  # xmin = peak.x - 500
  # xmax = peak.x + 500
  # ymin = peak.y - 500
  # ymax = peak.y + 500
  # bin_size = 50
  # xbin = int((xmax-xmin)/bin_size)
  # ybin = int((ymax-ymin)/bin_size)
  # h[f'xy_{ipeak+1}'] = ROOT.TH2D(f'xy_{ipeak+1}', f'xy_{ipeak+1}; x [#mum]; y [#mum]', xbin, xmin, xmax, ybin, ymin, ymax)

  startingPlate = peak.start
  if startingPlate < 5: continue
  endingPlate = peak.end
  h[f'theta_{ipeak+1}'] = ROOT.TH1D(f'theta_{ipeak+1}', f'theta_{ipeak+1}; theta [#mrad]',150, 0, 150)
  # for i in range(npl):
  #   idplate = ss.GetID(i)
  #   plate = idplate.ePlate
  #   nplate = ss.GetPlate(plate)
  for seg in segments:
    plate = seg.p
  # for plate in range(1,58):
    if plate < startingPlate or plate > endingPlate: continue
    x = seg.x
    y = seg.y
    if isInCluster(peak, x, y):
      tx = seg.tx
      ty = seg.ty
      vecs.append([tx, ty])
      theta =seg.theta*1000
      h[f'theta_{ipeak+1}'].Fill(theta)

    # ect = ROOT.EdbCouplesTree()
    # ect.InitCouplesTree("couples",datapath+f"/p{plate:03}/{brick}.{plate}.0.0.{cell}.cp.root","READ")
    # ect.eCut = cut
    # cutlist = ect.InitCutList()
    # if not ect.eTree: continue
    # if not cutlist:
    #   print("We have no entries, quitting!")
    #   continue
    # nsegcut = cutlist.GetN()
    # for ientry in range(nsegcut):
#       iseg = cutlist.GetEntry(ientry)
#       ect.GetEntry(iseg)
#       seg=ect.eS
# #      event = seg.MCEvt()
#       seg.SetZ(nplate.Z())
#       seg.SetPID(i)
#       seg.Transform(nplate.GetAffineXY())

        # sx = seg.X()
        # sy = seg.Y()
        # stx = seg.TX()
        # sty = seg.TY()
        # append_list_to_vecs(ipeak+1, [stx, sty])
        # ntuple.Fill(brick, ipeak+1, plate, sx, sy, stx, sty , stheta)
        # h[f'xy_{ipeak+1}'].Fill(sx,sy)
    # ect.Close()
  sph, direction = spherocity(vecs)
  mean = h[f'theta_{ipeak+1}'].GetMean()
  h['theta_mean'].Fill(mean)
  h['spherocity'].Fill(sph)
  h['direction'].Fill(direction[0])
  ntuple.Fill(cell, ix, iy, ipeak+1, peak.peak, peak.maxplate, startingPlate, endingPlate, peak.nseg, mean, sph, direction[0], peak.rankbin)

rootfile.cd()
# for event in vecs:
  # h[f'xy_{event}'].Write()
  # h[f'theta_{event}'].Write()
  # sph, direction = spherocity(vecs[event])
  # ntuple.Fill(cell,event,mean,sph)
h['theta_mean'].Write()
h['spherocity'].Write()
h['direction'].Write()
ntuple.Write()
rootfile.Close()
