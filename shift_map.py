import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max

ROOT.gStyle.SetOptStat("e")

def makeMap(hmap, tag):
    smoothed_hist = gaussian_filter(hmap, sigma=0.5, radius=1)
    peaks = peak_local_max(smoothed_hist, min_distance=2, threshold_abs=100)

    plt.figure(figsize=(8, 6))
    plt.imshow(smoothed_hist.T, origin='lower',cmap='viridis', interpolation='nearest')
    plt.plot(peaks[:, 0], peaks[:, 1], 'r.', markersize=10, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('2D Histogram with Detected Peaks')
    plt.xlabel('TX')
    plt.ylabel('TY')
    plt.legend()
    plt.savefig(f'/Users/fabioali/Desktop/shift_map_{tag}.png')

    return peaks


# Parameters
shiftRange = 50  #mrad
shiftStep  = 2   #mrad
nbins      = int((shiftRange * 2)/shiftStep) + 1 
nTag = 20

# File
path = '/Users/fabioali/cernbox/peaks_shift_1000.root'
file = ROOT.TFile.Open(path)
showers = file.Get("showers")

x = np.zeros(shiftRange+1)
y = np.zeros(shiftRange+1)

hmap = {}

for tag in range(1, nTag+1):
    hmap[tag], xedges, yedges = np.histogram2d(x, y, nbins, [[-shiftRange, shiftRange+2], [-shiftRange, shiftRange+2]])

for entry in showers:
    combination = int(entry.combination) - 1 
    peak = int(entry.rankbin)
    tag = showers.tag
    ix = combination // 31 
    iy = combination % 31 
    hmap[tag][ix, iy] = peak

for tag in range(1, nTag+1):
    peaks = makeMap(hmap[tag], tag)
    print(peaks)