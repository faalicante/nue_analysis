import ROOT
import os

outfile = ROOT.TFile.Open('signal_sample.root', 'RECREATE')
ntuple = ROOT.TNtuple('showers', 'Signal sample', 'event:energy:tx:ty:plate')

dz_shrink = 1315/1350
sampling = 0

with open ('/eos/experiment/sndlhc/users/falicant/shift_nue_Euniform_FLUKA_tuned/nue_int_10k.txt', 'r') as f:
    for line in f.readlines():
        line = line.strip().split(",")
        event = int(line[0])
        x0 = float(line[1])
        y0 = float(line[2])
        tx = float(line[5]) * dz_shrink
        ty = float(line[6]) * dz_shrink
        energy = float(line[7])
        plate = int(line[8]) - 3
        theta = ROOT.TMath.Sqrt(tx**2 + ty**2)
        if energy > 80 and plate < 40 and theta < 25 and theta > 5:
            print(f'./run.exe 1 {event} --x0 {x0} --y0 {y0} --tx {tx} --ty {ty} --p0 {plate}')
            os.system(f'./run.exe 1 {event} --x0 {x0} --y0 {y0} --tx {tx} --ty {ty} --p0 {plate}')
            ntuple.Fill(event, energy, tx, ty, plate)
            sampling += 1
        if sampling >= 10:
            break

outfile.cd()
ntuple.Write()
outfile.Write()
outfile.Close()