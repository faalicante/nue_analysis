import ROOT
from array import array

muon_file = "/eos/user/d/dannc/MuonBack_sim/03032022/sndLHC.Ntuple-TGeant4-1E5cm2.root"
muon_file_new = "/eos/experiment/sndlhc/users/dancc/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_z289.374023_BRICK11/9790422/sndLHC.Ntuple-muBkg1e5cm2_B11.root"
nue_file = "/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/sndLHC.Genie-TGeant4.root"
# file = ROOT.TFile.Open(muon_file)
file = ROOT.TFile.Open(nue_file)
cbmsim = file.cbmsim

# output_file = ROOT.TFile.Open("muon_showers.root","RECREATE")
output_file = ROOT.TFile.Open("nue_showers.root","RECREATE")
tree = ROOT.TTree("cbmsim","shower event")

N = 20100
_event = array('i', [0])
_count = array('i', [0])
_x = array('f', N*[0])
_y = array('f', N*[0])
_z = array('f', N*[0])
_tx = array('f', N*[0])
_ty = array('f', N*[0])
tree.Branch("event", _event, "event/I")
tree.Branch("count", _count, "count/I")
tree.Branch("x", _x, "x[count]/F")
tree.Branch("y", _y, "y[count]/F")
tree.Branch("z", _z, "z[count]/F")
tree.Branch("tx", _tx, "tx[count]/F")
tree.Branch("ty", _ty, "ty[count]/F")

for ievt in range(cbmsim.GetEntries()):
    count=0
    progress = ievt/cbmsim.GetEntries()*100
    if (int(progress%10))==0: print(int(progress), '%') 
    cbmsim.GetEntry(ievt)
    for track in cbmsim.MCTrack:
        x = track.GetStartX()
        y = track.GetStartY()
        z = track.GetStartZ()
        px = track.GetPx()
        py = track.GetPy()
        pz = track.GetPz()
        tx = px/pz
        ty = py/pz
        # if ROOT.TMath.Abs(z-307.8)>4.1: continue
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