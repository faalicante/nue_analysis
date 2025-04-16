import ROOT
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cell",dest="cell", type=int, required=True)
options = parser.parse_args()
cell = options.cell

ROOT.gStyle.SetOptStat("e")
map = {}

def makeTMap(hTmap, peaks):
    smoothed_hist = gaussian_filter(hTmap, sigma=1)
    # peaks = peak_local_max(smoothed_hist, min_distance=1, threshold_abs=100)
    peak_x_coords = peaks[0]
    peak_y_coords = peaks[1]
    peak_tx_coords = []
    peak_ty_coords = []
    output = ROOT.TFile.Open(f'peak_map_{cell}.root', 'RECREATE')
    outtree = ROOT.TNtuple("showers", "showers", "cell:combination:tx:ty:x:y")
    for px, py in zip(peak_x_coords, peak_y_coords):
        # print(px,py)
        x = np.searchsorted(xedges, px)
        y = np.searchsorted(yedges, py)
        combination = map[(x,y)][1]
        ix = combination // nbins
        iy = combination % nbins
        if abs(ix-25)==4 or abs(iy-25)==4:
            print(map[(x,y)], ix-25, iy-25, px, py)
        outtree.Fill(map[(x,y)][0], map[(x,y)][1], ix-25, iy-25, px, py)
        tx_coord = txedges[ix]
        ty_coord = tyedges[iy]
        peak_tx_coords.append(int(tx_coord)+1)
        peak_ty_coords.append(int(ty_coord)+1)
    plt.figure(figsize=(8, 6))
    plt.imshow(smoothed_hist.T, origin='lower',extent=[-shiftRange, shiftRange+shiftStep, -shiftRange, shiftRange+shiftStep],cmap='viridis', interpolation='none')
    plt.plot(peak_tx_coords, peak_ty_coords, 'r.', markersize=10, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('2D Histogram with Detected Peaks')
    plt.xlabel('TX')
    plt.ylabel('TY')
    # plt.legend()
    # plt.show()
    # plt.savefig(f'/Users/fabioali/Desktop/shifts/shift_Tmap_{cell}.png')
    plt.close()
    output.cd()
    outtree.Write()
    output.Close()

    return peak_tx_coords, peak_ty_coords

def makeMap(hmap):
    smoothed_hist = gaussian_filter(hmap, sigma=1)
    peaks = peak_local_max(smoothed_hist, min_distance=2, threshold_abs=0)
    peak_x_coords = xedges[peaks[:, 0]]
    peak_y_coords = yedges[peaks[:, 1]] 
    plt.figure(figsize=(8, 6))
    plt.imshow(smoothed_hist.T, origin='lower',extent=[xOffset,xRange+xOffset,yOffset,xRange+yOffset],cmap='viridis', interpolation='none')
    plt.plot(peak_x_coords+xStep/2, peak_y_coords+xStep/2, 'r.', markersize=5, label='Peaks')
    # plt.plot(xedg,yedg, 'r.', markersize=5, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('2D Histogram with Detected Peaks')
    plt.xlabel('X')
    plt.ylabel('Y')
    # plt.legend()
    # plt.show()
    # plt.savefig(f'/Users/fabioali/Desktop/shifts/shift_map_{cell}.png')
    plt.close()

    return peak_x_coords, peak_y_coords


# Parameters
shiftRange = 50  #mrad
shiftStep  = 2   #mrad
nbins      = int((shiftRange * 2)/shiftStep) + 1 
xRange = 12000  #um
xcell = cell%18
ycell = cell//18
xOffset = xcell*10000+4000       #um
yOffset = ycell*10000+4000      #um
# xOffset = 124000      #um
# yOffset = 64000      #um
xStep  = 150     #um
nxbins = int(xRange/xStep)
nTag = 20

# File
path = f'/eos/experiment/sndlhc/users/falicant/RUN1/b121/shift/{cell}/peaks.root'
# path = f'/Users/fabioali/cernbox/test_shift/trk/{cell}/peaks.root'
file = ROOT.TFile.Open(path)
showers = file.Get("showers")

tx = np.zeros(nbins)
ty = np.zeros(nbins)
x = np.zeros(nxbins)
y = np.zeros(nxbins)
ex = np.zeros(nxbins)
ey = np.zeros(nxbins)

# hTmap = {}

hmap, xedges, yedges = np.histogram2d(x, y, nxbins, range=[[xOffset,xRange+xOffset], [yOffset,xRange+yOffset]])
hTmap, txedges, tyedges = np.histogram2d(tx, ty, nbins, range=[[-shiftRange, shiftRange+shiftStep], [-shiftRange, shiftRange+shiftStep]])
# for tag in range(1, nTag+1):
#     hTmap[tag], txedges, tyedges = np.histogram2d(tx, ty, nbins, range=[[-shiftRange, shiftRange+shiftStep], [-shiftRange, shiftRange+shiftStep]])

for entry in showers:
    combination = int(entry.combination)
    if combination >2600: continue
    peak = int(entry.rankbin)
    tag = int(entry.tag)
    if tag > nTag : continue
    if peak<400: continue
    # if abs(entry.x-54000)>1000 or abs(entry.y-45000)>1000:continue
    x = (int(entry.x) - xOffset) // xStep
    y = (int(entry.y) - yOffset) // xStep
    ix = combination // nbins
    iy = combination % nbins
    # if abs(x-2)<3 or abs(x-nxbins+2)<3 or abs(y-2)<3 or abs(y-nxbins+2)<3 : continue
    # if ROOT.TMath.Sqrt((ix-25)**2+(iy-25)**2)<7: continue
    # if combination != 1096: continue
    if peak > hTmap[ix, iy]:
    # if combination == 702:
        hTmap[ix, iy] = peak
    if peak > hmap[x,y]:
        hmap[x,y] = peak
        map[(x, y)] = (entry.cell, combination)

# for tag in range(1, nTag+1):
    # peaksT = makeTMap(hTmap[tag], tag)
peaks = makeMap(hmap)
# print(len(peaks[0]))
peaksT = makeTMap(hTmap, peaks)
# print(len(peaksT[0]))
