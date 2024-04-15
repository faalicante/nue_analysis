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
    _sliceMax=c_int(0)
    _sliceMay=c_int(0)
    _sliceMaz=c_int(0)
    file[f'h2_{plate}'].GetXaxis().SetRangeUser(11500,14500)
    file[f'h2_{plate}'].GetYaxis().SetRangeUser(104000,107500)
    file[f'h2_{plate}'].GetMaximumBin(_sliceMax, _sliceMay, _sliceMaz)
    sliceMax=_sliceMax.value
    sliceMay=_sliceMay.value
    hx = file[f'h2_{plate}'].ProjectionX("",sliceMay-5,sliceMay+5,"")
    maxBinX = hx.GetMaximumBin()
    maxValueX = hx.GetBinCenter(maxBinX)
    fx = ROOT.TF1('fx','gaus')
    fitX = hx.Fit(fx,'SQ','',maxValueX-250,maxValueX+250)
    integralX = fx.Integral(maxValueX-250,maxValueX+250)

    fitMeanX = fitX.Parameter(1)
    fitVarX = fitX.Parameter(2)

    hy = file[f'h2_{plate}'].ProjectionY("",sliceMax-5,sliceMax+5,"")
    maxBinY = hy.GetMaximumBin()
    maxValueY = hy.GetBinCenter(maxBinY)
    fy = ROOT.TF1('fy','gaus')
    fitY = hx.Fit(fy,'SQ','',maxValueY-250,maxValueY+250)
    integralY = fy.Integral(maxValueY-250,maxValueY+250)
    fitMeanY = fitY.Parameter(1)
    fitVarY = fitY.Parameter(2)

    return fitMeanX, fitMeanY, fitVarX, fitVarY, integralX, integralY


couplesFile = 'hist_couples_aligned_1_11.root'
couples = loadHists(couplesFile)

z={57:0.00 ,56:-1340.35,55:-2714.28,54:-4036.24,53:-5337.83,52:-6684.41,51:-8001.78, 50:-9322.46 ,49:-10655.76, 48:-12021.93, 47:-13412.77, 46: -14795.63, 45:-16162.66, 44:-17534.37, 43:-18891.43, 42:-20246.43, 41:-21566.47 ,40:-22913.61 ,39:-24230.69 ,38:-25575.51 ,
   37:-26879.11 ,36:-28229.82 ,35:-29581.09 ,34:-30901.87 ,33:-32245.53 ,32:-33597.76 ,31:-34976.27 ,30:-36277.15 ,29:-37604.21 ,28:-38917.36 ,27:-40238.69, 26:-41584.73 ,25:-42908.92 ,24:-44276.78 ,23:-45643.12 ,22:-46972.20 ,21:-48355.71 ,20:-49722.71 ,19:-51066.87 ,
   18:-52427.69 ,17:-53749.28 ,16:-55136.12 ,15:-56493.45 ,14:-57828.26 ,13:-59199.22 ,12:-60574.84 ,11:-61938.31 ,10:-63250.14 , 9:-64578.07 , 8:-65892.94 , 7:-67241.09 , 6:-68555.77 , 5:-69898.01 , 4:-71225.37 , 3:-72564.90 , 2:-73891.68 , 1:-75247.02}
for key in range(1,25):
    z.pop(key, None)
for key in range(56,58):
    z.pop(key, None)
zList = sorted(z.values())
meanX = []
meanY = []
varX = []
varY = []

for plate in range(25, 56):
    fitMeanX, fitMeanY, fitVarX, fitVarY, integralX, integralY = fitCentroid(couples, plate)
    meanX.append(fitMeanX)
    meanY.append(fitMeanY)
    varX.append(fitVarX)
    varY.append(fitVarY)
    print(f'Centroid fit for plate {plate}: x {int(fitMeanX)} #pm {int(fitVarX)}, y {int(fitMeanY)} #pm {int(fitVarY)}, total width {int(integralX+integralY)}')

arrZ = array('d', zList)
arrX = array('d', meanX)
arrY = array('d', meanY)
errX = array('d', varX)
errY = array('d', varY)
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