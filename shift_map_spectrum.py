import ROOT

ROOT.gStyle.SetOptStat("e")

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cell",dest="cell", type=int, required=True)
options = parser.parse_args()
cell = options.cell


# Parameters
shiftRange = 50  #mrad
shiftStep  = 2   #mrad
nbins      = int((shiftRange * 2)/shiftStep) + 1 
xRange = 12000  #um
xcell = cell%18
ycell = cell//18
xOffset = xcell*10000+3000       #um
yOffset = ycell*10000+3000      #um
# xOffset = 44000       #um
# yOffset = 44000      #um
xStep  = 300     #um
nxbins = int(xRange/xStep)
nTag = 20
radius = 300
r0 = radius/xStep

def makeMap(hmap):
    spectrum = ROOT.TSpectrum2(2*nTag)
    # hmap.Smooth()
    # spectrum.Background(hmap, 10)
    nPeaks = spectrum.Search(hmap, 2, "noback")
    # nPeaks = spectrum.GetNPeaks()
    xPeaks = spectrum.GetPositionX()
    yPeaks = spectrum.GetPositionY()
    markers = []
    print(nPeaks)
    for i in range(nPeaks):
        x = xPeaks[i]
        y = yPeaks[i]
        # print(f"Peak at (x, y) = ({x}, {y})")
        marker = ROOT.TMarker(x, y, 20) 
        marker.SetMarkerColor(ROOT.kRed) 
        marker.SetMarkerSize(1) 
        marker.Draw()
        markers.append(marker)
    return markers


# def getMax(h2):
#     rankbin = h2.GetMaximum()
#     MaxBin = h2.GetMaximumBin()
#     ix,iy,iz
#     h2.GetBinXYZ(MaxBin, ix, iy, iz)
#     x = h2.GetXaxis().GetBinCenter(ix)
#     y = h2.GetYaxis().GetBinCenter(iy)
#     for(int iix = ix-r0 iix<=ix+r0 iix++):
#         for(int iiy = iy-r0 iiy<=iy+r0 iiy++):
#         h2.SetBinContent(iix,iiy,0)
#     return rankbin

# def get_peaks(h2, int *ranks):
#     TH2F *h2new = (TH2F*)h2.Clone("get_peaks")
#     for(int i=0 i<nTag i++):
#         int rankbin = getMax(*h2new, peaks, txt, bkg)
#         ranks[i] = rankbin


path = '/Users/fabioali/cernbox/test_shift'
# path = f'/eos/experiment/sndlhc/users/falicant/RUN1/b121/shift/{cell}/tag'
rootfile = ROOT.TFile(path+f"/shift_map.root","RECREATE")
ntuple = ROOT.TNtuple("shifts", "Tree of couples","x:y:combination:tx:ty:peak")
# File
peakfile = path+'/58/peaks.root'
file = ROOT.TFile.Open(peakfile)
showers = file.Get("showers")

hmap = ROOT.TH2D("hmap", "hmap;x;y", nxbins, xOffset, xOffset+xRange, nxbins, yOffset, yOffset+xRange)
hTmap = ROOT.TH2D("hTmap", "hTmap;TX;TY", nbins, -shiftRange, shiftRange, nbins, -shiftRange, shiftRange)
map = {}
for entry in showers:
    combination = int(entry.combination)
    # if entry.x<xOffset+2*xStep or entry.x>xOffset+xRange-2*xStep or entry.y<yOffset+2*xStep or entry.y>yOffset+xRange-2*xStep or (abs(entry.x-124750)<300 and abs(entry.y-52500)<1000): continue
    # if entry.x<xOffset+2*xStep or entry.x>xOffset+xRange-2*xStep or entry.y<yOffset+2*xStep or entry.y>yOffset+xRange-2*xStep: continue
    nseg = int(entry.nseg)
    rankbin = int(entry.rankbin)
    if rankbin<500:continue
    tag = int(entry.tag)
    if tag > nTag : continue
    x = entry.x
    y = entry.y
    Xbin = hmap.GetXaxis().FindBin(x)
    Ybin = hmap.GetYaxis().FindBin(y)
    # print(Xbin)
    peak = rankbin
    if peak > hmap.GetBinContent(Xbin, Ybin):
        hmap.SetBinContent(Xbin, Ybin, peak)
        map[(Xbin, Ybin)] = combination

c = ROOT.TCanvas("c","c",1100,600)
c.Divide(2,1)
c.cd(1)
xyPeaks = makeMap(hmap)
hmap.GetZaxis().SetRangeUser(500,hmap.GetMaximum())
hmap.Draw("colz")
# print(map)
for peak in xyPeaks:
    binPeak = hmap.FindBin(peak.GetX(), peak.GetY())
    Xbin = hmap.GetXaxis().FindBin(peak.GetX())
    Ybin = hmap.GetYaxis().FindBin(peak.GetY())
    if (Xbin, Ybin) not in map: continue
    combination = map[(Xbin, Ybin)]
    ix = combination // nbins + 1 
    iy = combination % nbins + 1
    tx = combination // nbins*2 - 50
    ty = combination % nbins*2 - 50
    # if ROOT.TMath.Sqrt((ix-26)**2+(iy-26)**2)<5: continue

    Tbin = hTmap.GetBin(ix,iy)

    nseg = hmap.GetBinContent(binPeak)
    ntuple.Fill(peak.GetX(), peak.GetY(), combination, tx, ty, nseg)
    if nseg > hTmap.GetBinContent(Tbin):
        hTmap.SetBinContent(Tbin, nseg)

c.cd(2)
hTmap.GetZaxis().SetRangeUser(500,hmap.GetMaximum())
hTmap.Draw("colz")
c.Draw()
c.Update()
rootfile.cd()
ntuple.Write()
c.Write()
rootfile.Write()
rootfile.Close()