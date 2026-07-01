import ROOT
import os
import re

def add_to_dat(event, xpos, ypos, projx, projy, ltx, lty, energy, plate):
  with open('nue_int_100.txt', 'a') as f:
    f.write(f'{event}, {xpos}, {ypos}, {projx}, {projy}, {ltx}, {lty}, {energy}, {plate}\n')


def decodeVolInt(wall_number, row_number, brick_number):
  wall = wall_number + 1
  column = brick_number
  row = row_number
  brick_map = {(0, 0): 2, (0, 1): 1, (1, 0): 4, (1, 1): 3}
  brick = brick_map[(row, column)]
  brick = wall*10 + brick
  return brick

def getVolInt(nu_vtx):
  nodeInt = ROOT.gGeoManager.FindNode(nu_vtx.X(), nu_vtx.Y(), nu_vtx.Z())
  pathInt = ROOT.gGeoManager.GetPath()

  wall_path = re.search(r"/Wall_(\d+)", pathInt)
  row_path = re.search(r"/Row_(\d+)", pathInt)
  brick_path = re.search(r"/Brick_(\d+)", pathInt)

  wall_number = int(wall_path.group(1)) if wall_path else None
  row_number = int(row_path.group(1)) if row_path else None
  brick_number = int(brick_path.group(1)) if brick_path else None
  if brick_number != None: brick = decodeVolInt(wall_number, row_number, brick_number)
  else: return None
  return brick

def getClosestEmuDetID(lep_track, emudetpoint):
  detID = -999
  if lep_track < -1: 
    print('Track not found')
    return detID
  for point in emudetpoint:
    if point.GetTrackID() == lep_track:
      detID = point.GetDetectorID()
      break
  return detID


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--from_evt", dest="from_evt", required=False, type=int, default=0)
parser.add_argument("--to_evt", dest="to_evt", required=True, type=int, default=None)
options = parser.parse_args()

from_evt = options.from_evt
to_evt = options.to_evt

path = '/eos/experiment/sndlhc/users/dancc/NUSIM/nue_inBrick21'

geoFile = path+'/1/geofile_full.Genie-TGeant4.root'
ROOT.TGeoManager.Import(geoFile)

# outputFile = ROOT.TFile('nue_int_100.root', 'RECREATE')
# outNtuple = ROOT.TNtuple('showers', 'showers', 'event:cell:xpos:ypos:ltx:lty:energy:plate')

sTree = ROOT.TChain('cbmsim')
for part in range(1,11):
  sTree.Add(path+f'/{part}/sndLHC.Genie-TGeant4.root')
# simName = 'inECC_sndLHC.Genie-TGeant4.root'
# simFile = ROOT.TFile.Open(simName)
# events = simFile.cbmsim

dz_shrink = 1315/1350

for i_event, event in enumerate(sTree):
  if i_event < from_evt: continue
  if i_event >= to_evt: break
  print(f'Processing event {i_event}')
  nutrack = event.MCTrack[0]
  leptrack = event.MCTrack[1]
  nu_vtx = ROOT.TVector3(nutrack.GetStartX(), nutrack.GetStartY(), nutrack.GetStartZ())
  lep_ang = ROOT.TVector3(leptrack.GetPx(), leptrack.GetPy(), leptrack.GetPz())
  nu_brick_int= getVolInt(nu_vtx)
  if nu_brick_int != 21: continue
  xem = nu_vtx.X() + 47.3 - 1.44
  yem = nu_vtx.Y() - 15.8 + 0.61
  xpos = round(xem * 1E4,2)
  ypos = round(yem * 1E4,2)
  ltx = round(lep_ang.X()/lep_ang.Z()*1000,2)
  lty = round(lep_ang.Y()/lep_ang.Z()*1000,2)
  energy = round(leptrack.GetEnergy(),2)
  detectorID = getClosestEmuDetID(1, event.EmulsionDetPoint)
  plate = detectorID % 1000
  projx = xpos - ltx/1000 * dz_shrink * ((plate-4) * 1350) ## -4 bc sim was made merging 60p neutrinos with 57p muons, skipping first 3 plates, i.e. p4 in neutrino is p1 in the mixed
  projy = ypos - lty/1000 * dz_shrink * ((plate-4) * 1350)
  add_to_dat(i_event, xpos, ypos, projx, projy, ltx, lty, energy, plate)
