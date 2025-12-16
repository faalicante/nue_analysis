import ROOT

file = ROOT.TFile.Open('peak_maps.root')
showers = file.showers

outfile = ROOT.TFile.Open('nue_tag.root', 'RECREATE')
ntuple = ROOT.TNtuple('ntuple', 'ntuple', 'event:dx:dy:tx:ty:dtx:dty:energy:plate:dplate:tag:peak:found')

with open ('nue_int_100.txt', 'r') as f:
    for line in f.readlines():
        line = line.strip().split(",")
        event = int(line[0])
        if event!=25:continue
        xProj = float(line[3])
        yProj = float(line[4])
        tx = float(line[5])
        ty = float(line[6])
        energy = float(line[7])
        plate = int(line[8])
        found = 0
        print(event)
        for shower in showers:
            if shower.cell != event: continue
            dx = shower.x - xProj
            dy = shower.y - yProj
            dist = ROOT.TMath.Sqrt(dx*dx + dy*dy)
            if dist > 200: continue
            dplate = shower.p - plate
            if dplate < 0: continue
            print("found")
            found = 1
            dtx = tx + shower.tx
            dty = ty + shower.ty
            tag = shower.tag
            peak = shower.peak
            ntuple.Fill(event, dx, dy, tx, ty, dtx, dty, energy, plate, dplate, tag, peak, found)
        if not found: 
            ntuple.Fill(event, -999, -999, tx, ty, -999, -999, energy, plate, -999, -999, -999, found)

ntuple.Write()
outfile.Write()
outfile.Close()