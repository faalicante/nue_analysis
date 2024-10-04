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

# path = '/eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b000022/hist_XYP_nue.root'
file = '/Users/fabioali/cernbox/hist_XYP_nue.root'
h = loadHists(file)
outputFile = ROOT.TFile("/Users/fabioali/Desktop/histo_combinations.root","RECREATE")

# SX(pl) = (1000+350)*pl*AX, SY(pl) = (1000+350)*pl*AY
# X(pl) = 1350*AX*pl + X0, Y(pl) = 1350*AY*pl + Y0
# (300/2)/(1350*56)=0.00198 150um shift = 3 bins
# +-30 mrad range, steps of 2mrad: 31x31 = 961 combinations

shiftSize  = 150 #um
binSize    = 50  #um
shiftRange = 30  #mrad
shiftStep  = 2   #mrad
stepZ      = 1000 + 350 #um

combination = 0
hComb = {}

xMin = 73000
xMax = 81000
yMin = 87000
yMax = 95000

xBin = int((xMax - xMin) / binSize)
yBin = int((yMax - yMin) / binSize)

c = ROOT.TCanvas("c", "c", 800, 800)
for shiftTX in np.arange(-shiftRange, shiftRange+1, shiftStep):
    for shiftTY in np.arange(-shiftRange, shiftRange+1, shiftStep):
        combination += 1
        if combination >1: continue
        print('Combination', combination)
        hComb[f'XY_{combination}'] = ROOT.TH2F(f'XY_{combination}', f'XY_{combination}', xBin, xMin, xMax, yBin, yMin, yMax)

        for layer in range(60):
            plate = layer + 1

            shiftX = shiftTX / 1000 * stepZ * layer
            shiftY = shiftTY / 1000 * stepZ * layer
            print(shiftTX, shiftTY, shiftX, shiftY)
            h[f'XYseg_{plate}'].GetXaxis().SetRangeUser(xMin, xMax)
            print('x', plate, h[f'XYseg_{plate}'].GetXaxis().GetXmin(), h[f'XYseg_{plate}'].GetXaxis().GetXmax())
            h[f'XYseg_{plate}'].GetYaxis().SetRangeUser(yMin, yMax)
            print('y', plate, h[f'XYseg_{plate}'].GetYaxis().GetXmin(), h[f'XYseg_{plate}'].GetYaxis().GetXmax())

            shiftedMinX = h[f'XYseg_{plate}'].GetXaxis().GetXmin() + shiftX
            shiftedMaxX = h[f'XYseg_{plate}'].GetXaxis().GetXmax() + shiftX
            h[f'XYseg_{plate}'].GetXaxis().SetLimits(shiftedMinX, shiftedMaxX)

            shiftedMinY = h[f'XYseg_{plate}'].GetYaxis().GetXmin() + shiftY
            shiftedMaxY = h[f'XYseg_{plate}'].GetYaxis().GetXmax() + shiftY
            h[f'XYseg_{plate}'].GetYaxis().SetLimits(shiftedMinY, shiftedMaxY)
            h[f'XYseg_{plate}'].Draw("colz")
            # c.Draw()
            # hComb[f'XY_{combination}'].Add(h[f'XYseg_{plate}']) ##check merge function
        print(h[f'XYseg_{plate}'].GetNbinsX(), hComb[f'XY_{combination}'].GetNbinsX())
        outputFile.cd()
        hComb[f'XY_{combination}'].Write()

outputFile.Close()