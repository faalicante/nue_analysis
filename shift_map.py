import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max


ROOT.gStyle.SetOptStat("e")

def makeTMap(hTmap, tag):
    smoothed_hist = gaussian_filter(hTmap, sigma=0.6)
    peaks = peak_local_max(hTmap, min_distance=2, threshold_abs=500)
    peak_x_coords = txedges[peaks[:, 0]] 
    peak_y_coords = tyedges[peaks[:, 1]] 
    plt.figure(figsize=(8, 6))
    plt.imshow(hTmap.T, origin='lower',extent=[-shiftRange, shiftRange+shiftStep, -shiftRange, shiftRange+shiftStep],cmap='viridis', interpolation='nearest')
    plt.plot(peak_x_coords+ shiftStep/2, peak_y_coords+ shiftStep/2, 'r.', markersize=10, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('2D Histogram with Detected Peaks')
    plt.xlabel('TX')
    plt.ylabel('TY')
    plt.legend()
    plt.savefig(f'/Users/fabioali/Desktop/shifts/shift_Tmap_{tag}.png')
    plt.close()

    return peak_x_coords, peak_y_coords

def makeMap(hmap, tag):
    smoothed_hist = gaussian_filter(hmap, sigma=1)
    peaks = peak_local_max(smoothed_hist, min_distance=5, threshold_abs=100)
    peak_x_coords = xedges[peaks[:, 0]]
    peak_y_coords = yedges[peaks[:, 1]] 
    plt.figure(figsize=(8, 6))
    plt.imshow(smoothed_hist.T, origin='lower',extent=[xOffset,xRange+xOffset,xOffset,xRange+xOffset],cmap='viridis', interpolation='nearest', vmin=0, vmax=200)
    plt.plot(peak_x_coords+xStep/2, peak_y_coords+xStep/2, 'r.', markersize=10, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('2D Histogram with Detected Peaks')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.savefig(f'/Users/fabioali/Desktop/shifts/shift_map_{tag}.png')
    plt.close()

    return peak_x_coords, peak_y_coords


# Parameters
shiftRange = 50  #mrad
shiftStep  = 2   #mrad
nbins      = int((shiftRange * 2)/shiftStep) + 1 
xRange = 200000  #um
xOffset = 0      #um
xStep  = 100     #um
nxbins      = int(xRange/xStep)
nTag = 20

# File
path = '/Users/fabioali/cernbox/peaks_shift.root'
file = ROOT.TFile.Open(path)
showers = file.Get("showers")

tx = np.zeros(nbins)
ty = np.zeros(nbins)
x = np.zeros(nxbins)
y = np.zeros(nxbins)

hTmap = {}

hmap, xedges, yedges = np.histogram2d(x, y, nxbins, range=[[xOffset,xRange+xOffset], [xOffset,xRange+xOffset]])
for tag in range(1, nTag+1):
    hTmap[tag], txedges, tyedges = np.histogram2d(tx, ty, nbins, range=[[-shiftRange, shiftRange+shiftStep], [-shiftRange, shiftRange+shiftStep]])

for entry in showers:
    combination = int(entry.combination) - 1 
    peak = int(entry.rankbin)
    tag = int(entry.tag)
    x = (int(entry.x) - xOffset) // xStep
    y = (int(entry.y) - xOffset) // xStep
    ix = combination // nbins
    iy = combination % nbins
    if peak > hTmap[tag][ix, iy]:
        hTmap[tag][ix, iy] = peak
    if peak > hmap[x,y]:
        hmap[x,y]= peak

for tag in range(1, nTag+1):
    peaksT = makeTMap(hTmap[tag], tag)
    peaks = makeMap(hmap[tag], tag)
    # print(len(peaks[0]))