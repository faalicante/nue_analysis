import ROOT
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max


ROOT.gStyle.SetOptStat("e")

def makeTMap(hTmap,):
    smoothed_hist = gaussian_filter(hTmap, sigma=1)
    peaks = peak_local_max(smoothed_hist, min_distance=2, threshold_abs=500)
    peak_x_coords = txedges[peaks[:, 0]] 
    peak_y_coords = tyedges[peaks[:, 1]] 
    plt.figure(figsize=(8, 6))
    plt.imshow(smoothed_hist.T, origin='lower',extent=[-shiftRange, shiftRange+shiftStep, -shiftRange, shiftRange+shiftStep],cmap='viridis', interpolation='none')
    # plt.plot(peak_x_coords+ shiftStep/2, peak_y_coords+ shiftStep/2, 'r.', markersize=10, label='Peaks')
    plt.colorbar(label='Counts')
    plt.title('Muon simulation')
    plt.xlabel('TX')
    plt.ylabel('TY')
    # plt.legend()
    plt.savefig(f'/Users/fabioali/Desktop/shifts/shift_Tmap_muon.png')
    plt.close()

    return peak_x_coords, peak_y_coords

def makeMap(hmap):
    smoothed_hist = gaussian_filter(hmap, sigma=1)
    peaks = peak_local_max(smoothed_hist, min_distance=5, threshold_abs=25)
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
    plt.savefig(f'/Users/fabioali/Desktop/shifts/shift_map_muon.png')
    plt.close()

    return peak_x_coords, peak_y_coords


# Parameters
shiftRange = 50  #mrad
shiftStep  = 2   #mrad
nbins      = int((shiftRange * 2)/shiftStep) + 1 
xRange = 13000  #um
xOffset = 288000      #um
yOffset = 84000      #um
xStep  = 100     #um
nxbins = int(xRange/xStep)
nTag = 10

# File
path = '/Users/fabioali/cernbox/peaks_shift_muon11.root'
file = ROOT.TFile.Open(path)
showers = file.Get("showers")
# pathSim = '/Users/fabioali/cernbox/nue_showers.root'
# fileSim = ROOT.TFile.Open(pathSim)
# cbmsim = fileSim.Get("cbmsim")

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
    combination = int(entry.combination) - 1 
    peak = int(entry.rankbin)
    tag = int(entry.tag)
    if tag > nTag : continue
    x = (int(entry.x) - xOffset) // xStep
    y = (int(entry.y) - yOffset) // xStep
    ix = combination // nbins
    iy = combination % nbins
    if peak > hTmap[ix, iy]:
        hTmap[ix, iy] = peak
    if peak > hmap[x,y]:
        hmap[x,y] = peak

# Example bin settings
# xedg = np.linspace(400, 0, 200000)
# yedg = np.linspace(400, 0, 200000)
# xedg = []
# yedg = []
# # "brick==22 && abs(ex+36.75)<9.25 && abs(ey-25.25)<9.25"
# for entry in cbmsim:
#     if entry.brick!=22:continue
#     ex = entry.ex
#     ey = entry.ey
#     if abs(ex+36.75)<9.25 and abs(ey-25.25)<9.25:
#         x = (ex+47.25) * 10000
#         y = (ey-15.75) * 10000
#         xedg.append(x)
#         yedg.append(y)
        

# for tag in range(1, nTag+1):
    # peaksT = makeTMap(hTmap[tag], tag)
peaksT = makeTMap(hTmap)
peaks = makeMap(hmap)
print(len(peaks[0]))
print(len(peaksT[0]))
# print(len(xedg))

#Brick_22: z=  307.7377cm  dZ=    3.8947cm  [  303.8430     311.6324] dx=    9.6498cm [  -45.9031     -26.6035] dy=    9.6490cm [   14.9472      34.2453]