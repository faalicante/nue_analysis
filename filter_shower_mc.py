import ROOT
from array import array
import progressbar

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
    file = ROOT.TFile.Open(nue_file)
    output_file = ROOT.TFile.Open(nue_path+"/nue_showers.root","RECREATE")
else:
    file = ROOT.TFile.Open(muon_file)
    output_file = ROOT.TFile.Open(nue_path+"/../muon1E5_simsndlhc/muon_showers.root","RECREATE")
cbmsim = file.cbmsim
tree = ROOT.TTree("cbmsim","shower event")

N = 100000
_brick = array('i', [0])
_event = array('i', [0])
_count = array('i', [0])
_nHits = array('i', [0])
_ex = array('f', [0])
_ey = array('f', [0])
_ez = array('f', [0])
_x = array('f', N*[0])
_y = array('f', N*[0])
_z = array('f', N*[0])
_tx = array('f', N*[0])
_ty = array('f', N*[0])
tree.Branch("brick", _brick, "brick/I")
tree.Branch("event", _event, "event/I")
tree.Branch("count", _count, "count/I")
tree.Branch("nHits", _nHits, "nHits/I")
tree.Branch("ex", _ex, "ex/F")
tree.Branch("ey", _ey, "ey/F")
tree.Branch("ez", _ey, "ey/F")
tree.Branch("x", _x, "x[count]/F")
tree.Branch("y", _y, "y[count]/F")
tree.Branch("z", _z, "z[count]/F")
tree.Branch("tx", _tx, "tx[count]/F")
tree.Branch("ty", _ty, "ty[count]/F")

bar = progressbar.ProgressBar(maxval=cbmsim.GetEntries(), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
bar.start()
for ievt in range(cbmsim.GetEntries()):
    bar.update(ievt+1)
    count=0
    brick_el = {11: 0, 12: 0, 13: 0, 14: 0, 21: 0, 22: 0, 23: 0, 24: 0, 31: 0, 32: 0, 33: 0, 34: 0, 41: 0, 42: 0, 43: 0, 44: 0, 51: 0, 52: 0, 53: 0, 54: 0, }
    cbmsim.GetEntry(ievt)

    for eHit in cbmsim.EmulsionDetPoint:
        trackID = eHit.GetTrackID()
        if trackID < 0: continue
        ex = cbmsim.MCTrack[1].GetStartX()
        ey = cbmsim.MCTrack[1].GetStartY()
        ez = cbmsim.MCTrack[1].GetStartZ()
        if cbmsim.MCTrack[trackID].GetProcID()==5 and cbmsim.MCTrack[trackID].GetEnergy() >= 0.1:
            detID = eHit.GetDetectorID()
            brickDet = DecodeBrickID(detID)
            brick_el[brickDet] += 1
            x = eHit.GetX()
            y = eHit.GetY()
            z = eHit.GetZ()
            px = eHit.GetPx()
            py = eHit.GetPy()
            pz = eHit.GetPz()
            tx = px/pz
            ty = py/pz
            _x[count] = x
            _y[count] = y
            _z[count] = z
            _tx[count] = tx
            _ty[count] = ty
            count+=1
    brickID = max(brick_el, key=brick_el.get)
    nHits = max(brick_el.values())
    if count>0:
        _brick[0] = brickID
        _event[0] = ievt
        _count[0] = count
        _nHits[0] = nHits
        _ex[0] = ex
        _ey[0] = ey
        _ez[0] = ez
        tree.Fill()
    # for track in cbmsim.MCTrack:
    #     x = track.GetStartX()
    #     y = track.GetStartY()
    #     z = track.GetStartZ()
    #     px = track.GetPx()
    #     py = track.GetPy()
    #     pz = track.GetPz()
    #     tx = px/pz
    #     ty = py/pz

    #     if track.GetProcID()==5 and track.GetEnergy() >= 0.1:
    #         _x[count] = x
    #         _y[count] = y
    #         _z[count] = z
    #         _tx[count] = tx
    #         _ty[count] = ty
    #         count+=1

bar.finish()
output_file.cd()
tree.Write()
output_file.Close()