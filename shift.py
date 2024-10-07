import ROOT
import numpy as np

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-p",dest="partition", type=int, required=True)
options = parser.parse_args()
partition = options.partition

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
            

# Nue simulation
path = '/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b000022'
# path = '/Users/fabioali/cernbox'
file = path + '/hist_XYP_nue.root'
h = loadHists(file)

combination = 0
# hComb = {}
combStart = partition * 100
combEnd = combStart + 100

for shiftTX in np.arange(-shiftRange, shiftRange+1, shiftStep):
    for shiftTY in np.arange(-shiftRange, shiftRange+1, shiftStep):
        combination += 1
        if combination < combStart or combination >= combEnd: continue
        # MC nue event in b22: 958, 445, 498, 1269 
        # if (shiftTX == 20 and shiftTY == -4) or (shiftTX == -26 and shiftTY == -10) or (shiftTX == 14 and shiftTY == 0) or (shiftTX == -4 and shiftTY == -16):
        outputFile = ROOT.TFile(path + f'/histo_shifts_{combination}.root',"RECREATE")
        print('Combination', combination)
        hComb = ROOT.TH2F('XYseg', 'XYseg', xBin, xMin, xMax, yBin, yMin, yMax)
        hList = ROOT.TList()
        for layer in range(60):
            plate = layer + 1
            shiftX = shiftTX / 1000 * stepZ * layer
            shiftY = shiftTY / 1000 * stepZ * layer
            hCrop = cropHist(h[f'XYseg_{plate}'], shiftX, shiftY)
            hCrop.Write()
            hComb.Add(hCrop)
        hComb.Write()
        outputFile.Close()