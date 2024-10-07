import ROOT
import numpy as np


def loadHists(histFile, query=None):
    f = ROOT.TFile.Open(histFile)
    histList = {}
    keyList = f.GetListOfKeys()
    for key in keyList:
        # if query is not None and key.GetName() not in query:
        #    continue
        hist = f.Get(key.GetName())
        hist.SetDirectory(ROOT.gROOT)
        hist.SetName(key.GetName())
        histList[key.GetName()] = hist
    if len(histList) == 0: raise Exception('ERROR: histList is empty!')
    f.Close()
    return histList


# SX(pl) = (1000+350)*pl*AX, SY(pl) = (1000+350)*pl*AY
# X(pl) = 1350*AX*pl + X0, Y(pl) = 1350*AY*pl + Y0
# (300/2)/(1350*56)=0.00198 150um shift = 3 bins
# +-30 mrad range, steps of 2mrad: 31x31 = 961 combinations

# Parameters
shiftSize  = 150 #um
binSize    = 50  #um
shiftRange = 30  #mrad
shiftStep  = 2   #mrad
stepZ      = 1000 + 350 #um

xMin = 71500
xMax = 82500
yMin = 85500
yMax = 96500

xBin = int((xMax - xMin) / binSize)
yBin = int((yMax - yMin) / binSize)

def cropHist(h2, shiftX, shiftY):
    xBin = int((xMax - xMin) / binSize)
    yBin = int((yMax - yMin) / binSize)
    hCrop = ROOT.TH2F(h2.GetTitle(), h2.GetTitle(), xBin, xMin, xMax, yBin, yMin, yMax)
    for xBin in range(1, h2.GetNbinsX()+1):
        xCenter = h2.GetXaxis().GetBinCenter(xBin) + shiftX
        if xCenter <= xMax and xCenter >= xMin:
            for yBin in range(1, h2.GetNbinsY()+1):
                yCenter = h2.GetYaxis().GetBinCenter(yBin) + shiftY
                if yCenter <= yMax and yCenter >= yMin:
                    content = h2.GetBinContent(xBin, yBin)
                    xBinNew = hCrop.GetXaxis().FindBin(xCenter)
                    yBinNew = hCrop.GetYaxis().FindBin(yCenter)
                    hCrop.SetBinContent(xBinNew, yBinNew, content)
    return hCrop
            

# path = '/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b000022/hist_XYP_nue.root'
file = '/Users/fabioali/cernbox/hist_XYP_nue.root'
h = loadHists(file)

combination = 0
hComb = {}

c = ROOT.TCanvas("c", "c", 800, 800)
for shiftTX in np.arange(-shiftRange, shiftRange+1, shiftStep):
    for shiftTY in np.arange(-shiftRange, shiftRange+1, shiftStep):
        combination += 1
        print('Combination', combination)
        hComb[f'XY_{combination}'] = ROOT.TH2F(f'XY_{combination}', f'XY_{combination}', xBin, xMin, xMax, yBin, yMin, yMax)
        hList = ROOT.TList()
        for layer in range(60):
            plate = layer + 1
            shiftX = shiftTX / 1000 * stepZ * layer
            shiftY = shiftTY / 1000 * stepZ * layer
            hCrop = cropHist(h[f'XYseg_{plate}'], shiftX, shiftY)
            hComb[f'XY_{combination}'].Add(hCrop)

# Saving histos
outputFile = ROOT.TFile("/Users/fabioali/Desktop/histo_combinations.root","RECREATE")
for combination in hComb.keys():
    print(combination)
    hComb[combination].Write()
outputFile.Close()