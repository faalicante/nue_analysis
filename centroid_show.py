import ROOT
import pandas as pd

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

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetTitleFontSize(0.1)
ROOT.gStyle.SetTitleSize(13,"XY")
ROOT.gStyle.SetTitleFont(43,"XY")
ROOT.gStyle.SetTitleOffset(1.15,"XY")
ROOT.gStyle.SetNdivisions(4,"XYZ")
ROOT.gStyle.SetLabelSize(0.04,"XYZ")

RUN = 1
WALL = 2
BRICK = 3

candidates = pd.read_csv(f"R{RUN}W{WALL}B{BRICK}_candidates.csv")
print(candidates)

for row in candidates.itertuples(index=False, name=None):
    cand = int(row[0])  
    xmin = row[3]
    xmax = row[4]
    ymin = row[5]
    ymax = row[6]
    xmin2 = (xmin + xmax)/2 - 2500
    xmax2 = (xmin + xmax)/2 + 2500
    ymin2 = (ymin + ymax)/2 - 2500
    ymax2 = (ymin + ymax)/2 + 2500
    

    # cellx = int(((xmin + xmax)/2) // 10000) +1 ##19 cells
    # celly = int(((ymin + ymax)/2) // 10000) +1
    cellx = int(((xmin + xmax)/2 + 5000) // 10000) ##18 cells
    celly = int(((ymin + ymax)/2 + 5000) // 10000)
    print(f"Candidate: {cand}, xmin: {xmin}, xmax: {xmax}, ymin: {ymin}, ymax: {ymax}, cellx: {cellx}, celly: {celly}")

    hist_file = f'/eos/user/f/falicant/nue_search/R{RUN}W{WALL}B{BRICK}/hist_couples_aligned_{cellx}_{celly}.root'
    hist_cp = loadHists(hist_file)

    c = ROOT.TCanvas("c","c",1920,1080)
    c.Divide(10,6)
    # i=0

    for p in range(1,58):
        # i+=1
        if p==16:continue
        c.cd(p)
        hist_cp[f'h_{p}'].SetTitle(f"Plate {p}")
        hist_cp[f'h_{p}'].GetYaxis().SetTitleOffset(1.25)
        hist_cp[f'h_{p}'].GetZaxis().SetRangeUser(0,40)
        hist_cp[f'h_{p}'].GetXaxis().SetRangeUser(xmin2,xmax2)
        hist_cp[f'h_{p}'].GetYaxis().SetRangeUser(ymin2,ymax2)
        hist_cp[f'h_{p}'].Draw("COLZ")
    c.SaveAs(f"R{RUN}W{WALL}B{BRICK}_{cand}.png")
'''

xmin = 70000
xmax = 71000
ymin = 83000
ymax = 84000
xmin2 = (xmin + xmax)/2 - 2500
xmax2 = (xmin + xmax)/2 + 2500
ymin2 = (ymin + ymax)/2 - 2500
ymax2 = (ymin + ymax)/2 + 2500

cellx = int(((xmin + xmax)/2) // 10000)
celly = int(((ymin + ymax)/2) // 10000)

hist_file = f'/eos/user/f/falicant/nue_search/R{RUN}W{WALL}B{BRICK}/hist_couples_aligned_{cellx}_{celly}.root'
print(hist_file)
print(xmin2,xmax2,ymin2,ymax2)
hist_cp = loadHists(hist_file)

c = ROOT.TCanvas("c","c",1920,1080)
c.Divide(3,2)
# i=0

for p in range(1,7):
    # i+=1
    c.cd(p)
    hist_cp[f'h_{p}'].SetTitle(f"Plate {p}")
    hist_cp[f'h_{p}'].GetYaxis().SetTitleOffset(1.25)
    hist_cp[f'h_{p}'].GetZaxis().SetRangeUser(0,30)
    hist_cp[f'h_{p}'].GetXaxis().SetRangeUser(xmin2,xmax2)
    hist_cp[f'h_{p}'].GetYaxis().SetRangeUser(ymin2,ymax2)
    hist_cp[f'h_{p}'].Draw("COLZ")
c.Draw()
'''