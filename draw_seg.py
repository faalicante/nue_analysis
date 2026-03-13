import ROOT
import ctypes

# from argparse import ArgumentParser
# parser = ArgumentParser()
# parser.add_argument('-e', '--event', type=int, required=True, help='Event to plot')
# args = parser.parse_args()

ROOT.gROOT.SetBatch(True)
stepZ = 1350
path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc/nuecc_muon_Euniform_RUN1/b000021'
opath = '/eos/user/f/falicant/shift/nue_regen'
# path = '/Users/fabioali/cernbox/shift/nue_regen'
# opath = path
ntuple = ROOT.TNtuple('showers', 'showers', 'event:comb:shiftTX:shiftTY:x:y:peak')
dz_shrink = 1315/1350

with open ('nue_int_100.txt', 'r') as f:
    for line in f.readlines():
        line = line.strip().split(",")
        event = int(line[0])
        # if event != args.event: continue
        xCenter = int(float(line[1]))
        yCenter = int(float(line[2]))
        tx = float(line[5]) * dz_shrink
        ty = float(line[6]) * dz_shrink
        binSize = 50
        xRange = 5000
        nbins = int(2*xRange/binSize)
        file = path+f'/b000021.0.0.{event+1}.trk.root'
        infile = ROOT.TFile.Open(file)
        tracks = infile.Get('tracks')
        c = ROOT.TCanvas('c', 'c', 800, 600)
        h = ROOT.TH2D('h', f'Event:{event}', nbins, xCenter-xRange, xCenter+xRange, nbins, yCenter-xRange, yCenter+xRange)
        shiftTX = -round(tx/2)*2
        shiftTY = -round(ty/2)*2
        if abs(shiftTX) > 50 or abs(shiftTY) > 50: combination = -999 
        else: combination = (shiftTY + 50)/2 * (50 + 1) + (shiftTX + 50)/2
        h_shift = []
        for layer in range(57):
            shiftX = shiftTX / 1000.0 * stepZ * layer
            shiftY = shiftTY / 1000.0 * stepZ * layer
            h_shift.append(ROOT.TH2D(f'h_shift_{layer}', f'Event:{event}', nbins, xCenter-xRange, xCenter+xRange, nbins, yCenter-xRange, yCenter+xRange))
            tracks.Draw(f's.eY+{shiftY}:s.eX+{shiftX}>>h_shift_{layer}', f's.eFlag==1&&s.eScanID.ePlate=={layer}')
            h.Add(h_shift[layer])
        h.Smooth()
        h.Draw('colz')
        c.SaveAs(f'{opath}/nue_shift_{event}.png')

        rankbin = h.GetMaximum() 
        maxBin = h.GetMaximumBin()
        ix = ctypes.c_int(0)
        iy = ctypes.c_int(0)
        iz = ctypes.c_int(0)
        h.GetBinXYZ(maxBin, ix, iy, iz)

        x = h.GetXaxis().GetBinCenter(ix.value)
        y = h.GetYaxis().GetBinCenter(iy.value)
        # print(f'Event: {event}, ShiftTX: {shiftTX}, ShiftTY: {shiftTY}, x: {x}, y: {y}, peak: {rankbin}')
        ntuple.Fill(event, combination, shiftTX, shiftTY, x, y, rankbin)

        del c
        del h
        infile.Close()

outfile = ROOT.TFile.Open(f'{opath}/nue_shifts.root', 'RECREATE')
ntuple.Write()
outfile.Write()
outfile.Close()