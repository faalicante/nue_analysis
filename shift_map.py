import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max


ROOT.gStyle.SetOptStat("e")

def makeTMap(hTmap):
    smoothed_hist = gaussian_filter(hTmap, sigma=0.6)
    peaks = peak_local_max(hTmap, min_distance=2, threshold_abs=150)
    peak_x_coords = txedges[peaks[:, 0]] 
    peak_y_coords = tyedges[peaks[:, 1]] 
    plt.figure(figsize=(8, 6))
    plt.imshow(smoothed_hist.T, origin='lower',extent=[-shiftRange, shiftRange+shiftStep, -shiftRange, shiftRange+shiftStep],cmap='viridis', interpolation='nearest')
    plt.plot(peak_x_coords+ shiftStep/2, peak_y_coords+ shiftStep/2, 'r.', markersize=10, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('2D Histogram with Detected Peaks')
    plt.xlabel('TX')
    plt.ylabel('TY')
    plt.legend()
    plt.savefig(f'/Users/fabioali/Desktop/shift_Tmap_test.png')

    return peak_x_coords, peak_y_coords

def makeMap(hmap):
    smoothed_hist = gaussian_filter(hmap, sigma=1)
    peaks = peak_local_max(smoothed_hist, min_distance=2, threshold_abs=150)
    peak_x_coords = xedges[peaks[:, 0]]
    peak_y_coords = yedges[peaks[:, 1]] 
    plt.figure(figsize=(8, 6))
    plt.imshow(smoothed_hist.T, origin='lower',extent=[72000,82000,86000,96000],cmap='viridis', interpolation='nearest')
    plt.plot(peak_x_coords+50, peak_y_coords+50, 'r.', markersize=10, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('2D Histogram with Detected Peaks')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.savefig(f'/Users/fabioali/Desktop/shift_map_test.png')

    return peak_x_coords, peak_y_coords


# Parameters
shiftRange = 30  #mrad
shiftStep  = 2   #mrad
nbins      = int((shiftRange * 2)/shiftStep) + 1 
nTag = 1

# File
path = '/Users/fabioali/cernbox/peaks_shift.root'
file = ROOT.TFile.Open(path)
showers = file.Get("showers")

tx = np.zeros(nbins)
ty = np.zeros(nbins)
x = np.zeros(100)
y = np.zeros(100)


hTmap, txedges, tyedges = np.histogram2d(tx, ty, nbins, range=[[-shiftRange, shiftRange+shiftStep], [-shiftRange, shiftRange+shiftStep]])
hmap, xedges, yedges = np.histogram2d(x, y, 100, range=[[72000,82000], [86000,96000]])

for entry in showers:
    combination = int(entry.comb) - 1 
    peak = int(entry.rankbin)
    x = (int(entry.x) - 72000) //100
    y = (int(entry.y) - 86000) //100
    ix = (combination // nbins) - shiftRange
    iy = (combination % nbins) - shiftRange
    contentT = hTmap[ix + shiftRange, iy + shiftRange]
    if peak > contentT:
        hTmap[ix + shiftRange, iy + shiftRange] = peak
    content = hmap[x,y]
    if peak > content:
        hmap[x,y] = peak

peaksT = makeTMap(hTmap)
peaks = makeMap(hmap)
print(peaksT)
print(peaks)