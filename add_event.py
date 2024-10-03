import ROOT
from array import array
import progressbar

def DecodeBrickID(detID):
    NWall = int(detID//1E4)
    NBrick = int((detID - NWall*1E4)//1E3)
    return int(f"{NWall}{NBrick}")

nue = 0 #1 to signal, 0 to muon bkg

muon_path = '/eos/user/d/dannc/MuonBack_sim/03032022'
nue_path = '/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022'
muon_file = muon_path+"/sndLHC.Ntuple-TGeant4-1E5cm2.root"
muon_file_new = "/eos/experiment/sndlhc/users/dancc/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_z289.374023_BRICK11/9790422/sndLHC.Ntuple-muBkg1e5cm2_B11.root"
nue_file = nue_path+"/sndLHC.Genie-TGeant4.root"

if nue:
    simfile = ROOT.TFile.Open(nue_file)
    peakfile = ROOT.TFile.Open(nue_path+'/sh_angle2.root')
    output_file = ROOT.TFile.Open(nue_path+"/nue_peaks.root","RECREATE")
else:
    simfile = ROOT.TFile.Open(muon_file)
    peakfile = ROOT.TFile.Open(nue_path+'/../muon1E5_simsndlhc/sh_angle2.root')
    output_file = ROOT.TFile.Open(nue_path+"/../muon1E5_simsndlhc/muon_peaks.root","RECREATE")
cbmsim = simfile.cbmsim
peaks = peakfile.couples
ntuple = ROOT.TNtuple("couples", "Tree of couples","brick:peak:max:maxpeak:start:end:nseg:theta:spherocity:dir:event:rank:energy")

for peak in peaks:
    print(f"Brick {peak.brick} Peak {peak.peak}")
    evt = int(peak.event)
    cbmsim.GetEntry(evt)
    energy = cbmsim.MCTrack[0].GetEnergy()
    ntuple.Fill(peak.brick, peak.peak, peak.max, peak.maxpeak, peak.start, peak.end, peak.nseg, peak.theta, peak.spherocity, peak.dir, peak.event, peak.rank, energy)

output_file.cd()
ntuple.Write()
output_file.Close()