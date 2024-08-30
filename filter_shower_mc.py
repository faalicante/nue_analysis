import ROOT
from array import array

def DecodeBrickID(detID):
    NWall = int(detID//1E4)
    NBrick = int((detID - NWall*1E4)//1E3)
    return int(f"{NWall}{NBrick}")

nue = 1 #1 to signal, 0 to muon bkg

muon_path = '/eos/user/d/dannc/MuonBack_sim/03032022'
nue_path = '/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022'
muon_file = muon_path+"/sndLHC.Ntuple-TGeant4-1E5cm2.root"
muon_file_new = "/eos/experiment/sndlhc/users/dancc/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_z289.374023_BRICK11/9790422/sndLHC.Ntuple-muBkg1e5cm2_B11.root"
nue_file = nue_path+"/sndLHC.Genie-TGeant4.root"

if nue:
    output_file = ROOT.TFile.Open(nue_path+"/nue_showers.root","RECREATE")
    file = ROOT.TFile.Open(nue_file)
else:
    output_file = ROOT.TFile.Open(muon_path+"/muon_showers.root","RECREATE")
    file = ROOT.TFile.Open(muon_file)
cbmsim = file.cbmsim
tree = ROOT.TTree("cbmsim","shower event")

N = 20100
_brick = array('i', [0])
_event = array('i', [0])
_count = array('i', [0])
_x = array('f', N*[0])
_y = array('f', N*[0])
_z = array('f', N*[0])
_tx = array('f', N*[0])
_ty = array('f', N*[0])
tree.Branch("brick", _brick, "brick/I")
tree.Branch("event", _event, "event/I")
tree.Branch("count", _count, "count/I")
tree.Branch("x", _x, "x[count]/F")
tree.Branch("y", _y, "y[count]/F")
tree.Branch("z", _z, "z[count]/F")
tree.Branch("tx", _tx, "tx[count]/F")
tree.Branch("ty", _ty, "ty[count]/F")

for ievt in range(cbmsim.GetEntries()):
    brick_el = {11: 0, 12: 0, 13: 0, 14: 0, 21: 0, 22: 0, 23: 0, 24: 0, 31: 0, 32: 0, 33: 0, 34: 0, 41: 0, 42: 0, 43: 0, 44: 0, 51: 0, 52: 0, 53: 0, 54: 0, }
    count=0
    progress = ievt/cbmsim.GetEntries()*100
    if (int(progress%10))==0: print(int(progress), '%') 
    cbmsim.GetEntry(ievt)

    for eHit in cbmsim.EmulsionDetPoint:
        detID = eHit.GetDetectorID()
        brickDecoded = DecodeBrickID(detID)
        brick_el[brickDecoded] += 1
        brickID = max(brick_el, key=brick_el.get)
    _brick[0] = brickID

    for track in cbmsim.MCTrack:
        x = track.GetStartX()
        y = track.GetStartY()
        z = track.GetStartZ()
        px = track.GetPx()
        py = track.GetPy()
        pz = track.GetPz()
        tx = px/pz
        ty = py/pz

        if track.GetProcID()==5 and track.GetEnergy() >= 0.1:
            _x[count] = x
            _y[count] = y
            _z[count] = z
            _tx[count] = tx
            _ty[count] = ty
            count+=1

    if count>0:
        _event[0] = ievt
        _count[0] = count
        tree.Fill()

output_file.cd()
tree.Write()
output_file.Close()