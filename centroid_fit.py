import ROOT
from array import array
from ctypes import c_int


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

def fitCentroid(file, plate):
    _binMaxX=c_int(0)
    _binMaxY=c_int(0)
    _binMaxZ=c_int(0)
    file[f'h_{plate}'].GetXaxis().SetRangeUser(11500+plate*10,13500+plate*10)
    file[f'h_{plate}'].GetYaxis().SetRangeUser(104000+plate*10,106000+plate*10)
    file[f'h_{plate}'].GetMaximumBin(_binMaxX, _binMaxY, _binMaxZ)
    binMaxX=_binMaxX.value
    binMaxY=_binMaxY.value
    
    hx = file[f'h_{plate}'].ProjectionX("",binMaxY-10,binMaxY+10,"")
    maxBinX = hx.GetMaximumBin()
    maxValueX = hx.GetBinCenter(maxBinX)
    fx = ROOT.TF1('fx','gaus(0)+pol0(3)')
    fx.SetParameter(1,maxValueX)
    fx.SetParameter(2,250)
    fitX = hx.Fit(fx,'SQ')
    integralX = fx.Integral(maxValueX-500,maxValueX+500)
    fitMeanX = fitX.Parameter(1)
    fitMeanErrX = fitX.ParError(1)
    fitVarX = fitX.Parameter(2)

    hy = file[f'h_{plate}'].ProjectionY("",binMaxX-10,binMaxX+10,"")
    maxBinY = hy.GetMaximumBin()
    maxValueY = hy.GetBinCenter(maxBinY)
    fy = ROOT.TF1('fy','gaus(0)+pol0(3)')
    fy.SetParameter(1,maxValueY)
    fy.SetParameter(2,250)
    fitY = hy.Fit(fy,'SQ')
    integralY = fy.Integral(maxValueY-500,maxValueY+500)
    fitMeanY = fitY.Parameter(1)
    fitMeanErrY = fitY.ParError(1)
    fitVarY = fitY.Parameter(2)

    return fitMeanX, fitMeanY, fitMeanErrX, fitMeanErrY, fitVarX, fitVarY, integralX, integralY


couplesFile = '/afs/cern.ch/work/f/falicant/public/nue_search/hist_couples_aligned_mos.root'
couples = loadHists(couplesFile)

#set for b21 cell_1_11
#z={57:0.00 ,56:-1340.35,55:-2714.28,54:-4036.24,53:-5337.83,52:-6684.41,51:-8001.78, 50:-9322.46 ,49:-10655.76, 48:-12021.93, 47:-13412.77, 46: -14795.63, 45:-16162.66, 44:-17534.37, 43:-18891.43, 42:-20246.43, 41:-21566.47 ,40:-22913.61 ,39:-24230.69 ,38:-25575.51 ,
#   37:-26879.11 ,36:-28229.82 ,35:-29581.09 ,34:-30901.87 ,33:-32245.53 ,32:-33597.76 ,31:-34976.27 ,30:-36277.15 ,29:-37604.21 ,28:-38917.36 ,27:-40238.69, 26:-41584.73 ,25:-42908.92 ,24:-44276.78 ,23:-45643.12 ,22:-46972.20 ,21:-48355.71 ,20:-49722.71 ,19:-51066.87 ,
#   18:-52427.69 ,17:-53749.28 ,16:-55136.12 ,15:-56493.45 ,14:-57828.26 ,13:-59199.22 ,12:-60574.84 ,11:-61938.31 ,10:-63250.14 , 9:-64578.07 , 8:-65892.94 , 7:-67241.09 , 6:-68555.77 , 5:-69898.01 , 4:-71225.37 , 3:-72564.90 , 2:-73891.68 , 1:-75247.02}

z = {57: 34655.96, 56: 33329.90, 55: 32006.92, 54: 30685.82, 53: 29372.04, 52: 28061.45, 51: 26715.19, 50: 25373.10, 49: 24062.98, 48: 22730.45, 47: 21352.18, 46: 19961.56,
    45: 18634.26, 44: 17268.02, 43: 15889.72, 42: 14523.08, 41: 13210.54, 40: 11906.15, 39: 10572.31, 38:  9258.64, 37:  7942.98, 36:  6613.56, 35:  5267.96, 34:  3965.65,
    33:  2663.28, 32:  1344.32, 31:0.00, 30: -1317.17, 29: -2660.18, 28: -3936.19, 27: -5261.42, 26: -6600.40, 25: -7903.31, 24: -9317.59, 23:-10680.85, 22:-11971.73,
    21:-13345.02, 20:-14675.77, 19:-16021.16, 18:-17358.80, 17:-18693.34, 16:-20052.59, 15:-21400.94, 14:-22716.23, 13:-24059.38, 12:-25422.78, 11:-26722.78, 10:-28046.87,
    9:-29377.34,  8:-30695.57,  7:-32045.26,  6:-33379.57,  5:-34737.75,  4:-36102.56,  3:-37462.51,  2:-38769.36,  1:-40109.79}
z.update((x, y-34655.96) for x, y in z.items())
for key in range(1,18):
     z.pop(key, None)
for key in range(52,58):
     z.pop(key, None)
zList = sorted(z.values())
meanX = []
meanY = []
meanErrX = []
meanErrY = []
varX = []
varY = []

for plate in range(18, 52):
    fitMeanX, fitMeanY, fitMeanErrX, fitMeanErrY, fitVarX, fitVarY, integralX, integralY = fitCentroid(couples, plate)
    meanX.append(fitMeanX)
    meanY.append(fitMeanY)
    meanErrX.append(fitMeanErrX)
    meanErrY.append(fitMeanErrY)
    varX.append(fitVarX)
    varY.append(fitVarY)
    print(f'Centroid fit for plate {plate}: x {int(fitMeanX)} #pm {int(fitVarX)}, y {int(fitMeanY)} #pm {int(fitVarY)}, total width {int(integralX+integralY)}')

arrZ = array('d', zList)
arrX = array('d', meanX)
arrY = array('d', meanY)
errX = array('d', meanErrX)
errY = array('d', meanErrY)
errZ = array('d', [0]*len(zList))
gx = ROOT.TGraphErrors(len(zList), arrZ, arrX, errZ, errX)
gy = ROOT.TGraphErrors(len(zList), arrZ, arrY, errZ, errY)
gx.SetTitle('Centroid X')
gx.GetXaxis().SetTitle('z[#mum]')
gx.GetYaxis().SetTitle('x[#mum]')
gx.GetYaxis().SetTitleOffset(1.4)
gy.SetTitle('Centroid Y')
gy.GetXaxis().SetTitle('z[#mum]')
gy.GetYaxis().SetTitle('y[#mum]')
gy.GetYaxis().SetTitleOffset(1.4)
fx = ROOT.TF1("fz","pol1")
fy = ROOT.TF1("fy","pol1")
gx.Fit(fx,"") 
gy.Fit(fy,"") 
slopeX = fx.GetParameter(1)
slopeErrX = fx.GetParError(1)
slopeY = fy.GetParameter(1)
slopeErrY = fy.GetParError(1)

c = ROOT.TCanvas('c','centroids', 1500, 600)
c.Divide(2,1)
c.cd(1)  
gx.SetMarkerColor(4)
gx.SetMarkerStyle(21)
gx.Draw("AP")
c.cd(2)
gy.SetMarkerColor(4)
gy.SetMarkerStyle(21)
gy.Draw("AP")