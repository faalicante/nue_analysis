import ROOT

file1 = ROOT.TFile.Open('nue_shifts.root')
showers1 = file1.showers

outfile = ROOT.TFile.Open('nue_tag5.root', 'RECREATE')
ntuple1 = ROOT.TNtuple('showers1', 'Signal only compared to MC', 'event:dx:dy:tx:ty:dtx:dty:energy:plate:found')
ntuple2 = ROOT.TNtuple('showers2', 'Full tag compared with signal only', 'event:dx:dy:tx:ty:dtx:dty:energy:plate:dplate:tag:peak:found:bkg:nseg')

dz_shrink = 1315/1350

with open ('nue_int_10k.txt', 'r') as f:
    for line in f.readlines():
        line = line.strip().split(",")
        event = int(line[0])
        # if event > 999: break
        if event%100 == 0: print(f'Event: {event}') 
        xProj = float(line[3])
        yProj = float(line[4])
        tx = float(line[5]) * dz_shrink
        ty = float(line[6]) * dz_shrink
        energy = float(line[7])
        plate = int(line[8]) - 3
        found = 0
        mindist = 1000
        # print(event)
        for shower1 in showers1:
            if shower1.event != event: continue
            if shower1.peak ==0:
                print(f'Event: {event} has no peak. Skipping.')
                found = -1
                break
            x1 = shower1.x
            y1 = shower1.y
            comb1 = shower1.comb
            dx1 = x1 - xProj
            dy1 = y1 - yProj
            dtx1 = tx + shower1.shiftTX
            dty1 = ty + shower1.shiftTY
        file2 = ROOT.TFile.Open(f'maps/peak_map_{event}.root')
        showers2 = file2.showers
        for shower2 in showers2:
            if shower2.cell != event: continue
            dx2 = shower2.x - x1
            dy2 = shower2.y - y1
            dist = ROOT.TMath.Sqrt(dx2*dx2 + dy2*dy2)
            if dist > 200: continue
            dplate = shower2.p - plate
            # if dplate < -5: continue
            dtx2 = tx + shower2.tx
            dty2 = ty + shower2.ty
            # if abs(dtx) > 10 or abs(dty) > 10: continue
            found = 1
            if dist < mindist:
                mindist = dist
                min_dx = dx2
                min_dy = dy2
                min_dtx = dtx2
                min_dty = dty2
                min_dplate = dplate
                tag = shower2.tag
                peak = shower2.peak
                nseg = shower2.nseg
                bkg = shower2.bkg
        ntuple1.Fill(event, dx1, dy1, tx, ty, dtx1, dty1, energy, plate, found)
        if found: ntuple2.Fill(event, min_dx, min_dy, tx, ty, min_dtx, min_dty, energy, plate, min_dplate, tag, peak, found, bkg, nseg)
        else: ntuple2.Fill(event, -999, -999, tx, ty, -999, -999, energy, plate, -999, -999, found, -999, -999)
        file2.Close()   

outfile.cd()
ntuple1.Write()
ntuple2.Write()
outfile.Write()
outfile.Close()
file1.Close()