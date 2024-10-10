import ROOT

ROOT.gStyle.SetOptStat("e")

# Parameters
shiftRange = 30  #mrad
shiftStep  = 2   #mrad
nbins      = int((shiftRange * 2)/shiftStep) + 1 

# File
path = '/Users/fabioali/cernbox/peaks_shift.root'
file = ROOT.TFile.Open(path)
showers = file.Get("showers")

hmap = ROOT.TH2F("hmap","hmap;TX;TY", nbins, -shiftRange, shiftRange+2, nbins, -shiftRange, shiftRange+2)

for entry in showers:
    combination = int(entry.comb)
    peak = int(entry.rankbin)
    ix = combination//31 +1
    iy = combination%31 +1 
    content = hmap.GetBinContent(ix,iy)
    if peak > content:
        hmap.SetBinContent(ix, iy, peak)

c = ROOT.TCanvas("c","c",800,800)
hmap.Draw("colz")
c.Print('/Users/fabioali/Desktop/shift_map.png')