import ROOT
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cell",dest="cell", type=int, required=True)
parser.add_argument("-d",dest="data", type=int, required=True)
options = parser.parse_args()
cell = options.cell
data = options.data

ROOT.gStyle.SetOptStat("e")
map = {}

def makeTMap(hTmap, peaks):
    smoothed_hist = gaussian_filter(hTmap, sigma=0)
    # peaks = peak_local_max(smoothed_hist, min_distance=1, threshold_abs=100)
    peak_x_coords = peaks[0]
    peak_y_coords = peaks[1]
    peak_tx_coords = []
    peak_ty_coords = []
    output = ROOT.TFile.Open(f'peak_map_{cell}.root', 'RECREATE')
    outtree = ROOT.TNtuple("showers", "showers", "cell:combination:tag:tx:ty:x:y:p:peak")
    for px, py in zip(peak_x_coords, peak_y_coords):
        x = np.searchsorted(xedges, px)
        y = np.searchsorted(yedges, py)
        if (x,y) not in map: continue
        combination = map[(x,y)][1]
        ix = combination % nbins
        iy = combination // nbins
        if abs(px-xCenter)>5000 or abs(py-yCenter)>5000: continue
        outtree.Fill(map[(x,y)][0], map[(x,y)][1], map[(x,y)][2], ix*2-50, iy*2-50, map[(x,y)][5], map[(x,y)][6], map[(x,y)][3], map[(x,y)][4])
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
    plt.legend()
    # plt.show()
    plt.savefig(f'shift_Tmap_{cell}.png')
    plt.close()
    output.cd()
    outtree.Write()
    output.Close()

    return peak_tx_coords, peak_ty_coords

def makeMap(hmap):
    smoothed_hist = gaussian_filter(hmap, sigma=0)
    peaks = peak_local_max(smoothed_hist, min_distance=2, threshold_abs=0)
    peak_x_coords = xedges[peaks[:, 0]]
    peak_y_coords = yedges[peaks[:, 1]]
    plt.figure(figsize=(8, 6))
    plt.imshow(smoothed_hist.T, origin='lower',extent=[xCenter-xRange, xCenter+xRange, yCenter-xRange, yCenter+xRange],cmap='viridis', interpolation='none')
    plt.plot(peak_x_coords+xStep/2, peak_y_coords+xStep/2, 'r.', markersize=2, label='Peaks')
    # plt.plot(xedg,yedg, 'r.', markersize=5, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('2D Histogram with Detected Peaks')
    plt.xlabel('X')
    plt.ylabel('Y')
    # plt.legend()
    # plt.show()
    plt.savefig(f'shift_map_{cell}.png')
    plt.close()

    return peak_x_coords, peak_y_coords


# Parameters
shiftRange = 50  #mrad
shiftStep  = 2   #mrad
nbins      = int((shiftRange * 2)/shiftStep) + 1
dz = 1350
if data==1:
    xRange = 6000
    with open ('nue_int_100.txt', 'r') as f:
        lines = f.readlines()
        line = lines[cell].strip().split(",")
        xCenter = float(line[3])
        yCenter = float(line[4])
else:
    xRange = 10000  #um
    xcell = cell%18
    ycell = cell//18
    xCenter = (xcell+1)*10000
    yCenter = (ycell+1)*10000
# xOffset = 288000      #um
# yOffset = 83000      #um
xStep  = 100     #um
nxbins = int(2*xRange/xStep)
nTag = 50

# File
# path = f'/eos/experiment/sndlhc/users/falicant/RUN1/b121/shift/{cell}/peaks.root'
# path = f'/Users/fabioali/cernbox/test_shift/trk/{cell}/peaks.root'
path = f'peaks_{cell}.root'
file = ROOT.TFile.Open(path)
showers = file.Get("showers")

tx = np.zeros(nbins)
ty = np.zeros(nbins)
x = np.zeros(nxbins)
y = np.zeros(nxbins)
ex = np.zeros(nxbins)
ey = np.zeros(nxbins)


hmap, xedges, yedges = np.histogram2d(x, y, nxbins, range=[[xCenter-xRange, xCenter+xRange], [yCenter-xRange, yCenter+xRange]])
hTmap, txedges, tyedges = np.histogram2d(tx, ty, nbins, range=[[-shiftRange, shiftRange+shiftStep], [-shiftRange, shiftRange+shiftStep]])

for entry in showers:
    if entry.nfound > 1: continue
    combination = int(entry.combination)
    peak = int(entry.rankbin)
    tag = int(entry.tag)
    if tag > nTag : continue
    # if peak<500: continue
    ix = combination % nbins
    iy = combination // nbins
    tx = ix*2-50
    ty = iy*2-50
    plate = entry.maxplate
    projx = entry.x - tx/1000 * ((plate-1) * dz)
    projy = entry.y - ty/1000 * ((plate-1) * dz)
    if abs(projx-xCenter)>xRange or abs(projy-yCenter)>xRange: continue
    x = int((projx - xCenter + xRange) // xStep)
    y = int((projy - yCenter + xRange) // xStep)
    # if (combination==1293 and tag==1):
    #     print(plate, tag, x, y, tx, ty, projx, entry.x, entry.y, peak, hmap[x,y], map[(x, y)])
    # if abs(x-2)<3 or abs(x-nxbins+2)<3 or abs(y-2)<3 or abs(y-nxbins+2)<3 : continue
    # if ROOT.TMath.Sqrt((ix-25)**2+(iy-25)**2)<7: continue
    # if combination != 1234: continue
    if peak > hTmap[ix, iy]:
        hTmap[ix, iy] = peak
    if peak > hmap[x,y]:
        hmap[x,y] = peak
        map[(x, y)] = (entry.cell, combination, tag, entry.start, peak, entry.x, entry.y)
# print(map[(67,59)])


peaks = makeMap(hmap)
peaksT = makeTMap(hTmap, peaks)