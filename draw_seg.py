import ROOT
import ctypes
import os

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-e', '--event', dest='event', type=int, required=False, help='Event to plot')
parser.add_argument('-p', '--partition', dest='partition', type=int, required=False, help='Partition of 100 events')
args = parser.parse_args()

part = args.partition
ROOT.gROOT.SetBatch(True)
stepZ = 1350
path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc/nuecc_muon_Euniform_RUN1_FLUKA25/b000021'
opath = '/eos/user/f/falicant/shift/nue_Euniform_FLUKA'
# path = '/Users/fabioali/cernbox/shift/nue_regen'
# opath = path

outfile = ROOT.TFile.Open(f'{opath}/nue_shifts_{part}.root', 'RECREATE')
ntuple = ROOT.TNtuple('showers', 'showers', 'event:comb:shiftTX:shiftTY:x:y:peak')
dz_shrink = 1315/1350

with open ('/eos/experiment/sndlhc/users/falicant/shift_nue_regen_100/nue_int_10k.txt', 'r') as f:
    for line in f.readlines():
        line = line.strip().split(",")
        event = int(line[0])
        if event < part*100:continue
        if event >= (part+1)*100:break
        # if event != args.event: continue
        xCenter = int(float(line[1]))
        yCenter = int(float(line[2]))
        tx = float(line[5]) * dz_shrink
        ty = float(line[6]) * dz_shrink
        binSize = 50
        xRange = 5000
        nbins = int(2*xRange/binSize)
        file = path+f'/b000021.0.0.{event+1}.trk.root'
        if not os.path.exists(file):
            print(f'File {file} does not exist. Skipping event {event}.')
            continue
        infile = ROOT.TFile.Open(file)
        tracks = infile.Get('tracks')
        if tracks.GetEntries("s.eFlag==1") == 0:
            print(f'Event {event} has no valid tracks. Skipping.')
            continue
        else:
            print(f'Event {event} has {tracks.GetEntries("s.eFlag==1")} valid tracks. Processing.')
        # h = ROOT.TH2D('h', f'Event:{event}', nbins, xCenter-xRange, xCenter+xRange, nbins, yCenter-xRange, yCenter+xRange)
        h_shift = ROOT.TH2D('h_shift', f'Event:{event}', nbins, xCenter-xRange, xCenter+xRange, nbins, yCenter-xRange, yCenter+xRange)
        shiftTX = -round(tx/2)*2
        shiftTY = -round(ty/2)*2
        if abs(shiftTX) > 50 or abs(shiftTY) > 50: combination = -999 
        else: combination = (shiftTY + 50)/2 * (50 + 1) + (shiftTX + 50)/2
        # tracks.Draw(f's.eY:s.eX>>h', f's.eFlag==1')
        h_shifts = []
        for layer in range(57):
            plate = layer + 1
            shiftX = shiftTX / 1000.0 * stepZ * layer
            shiftY = shiftTY / 1000.0 * stepZ * layer
            h_shifts.append(ROOT.TH2D(f'h_shift_{plate}', f'Event:{event}', nbins, xCenter-xRange, xCenter+xRange, nbins, yCenter-xRange, yCenter+xRange))
            tracks.Draw(f's.eY+{shiftY}:s.eX+{shiftX}>>h_shift_{plate}', f's.eFlag==1&&s.eScanID.ePlate=={plate}')
            h_shift.Add(h_shifts[layer])
        # c = ROOT.TCanvas('c', 'c', 1500, 600)
        # c.Divide(2,1)
        # h_shift.Smooth()
        # c.cd(1)
        # h.Draw('colz')
        # c.cd(2)
        # h_shift.Draw('colz')
        # c.SaveAs(f'{opath}/nue_shift2_{event}.png')

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

ntuple.Write()
outfile.Write()
outfile.Close()