import ROOT
import subprocess
import re

outfile = ROOT.TFile.Open('signal_sample.root', 'RECREATE')
ntuple = ROOT.TNtuple('showers', 'Signal sample', 'event:energy:tx:ty:plate:zScale')

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
            print(f'./run.exe 1 {event} --x0 {x0} --y0 {y0} --tx {tx:.2f} --ty {ty:.2f} --p0 {plate}')


            cmd = [
                "./run.exe",
                "1",
                str(event),
                "--x0", str(x0),
                "--y0", str(y0),
                "--tx", f"{tx:.2f}",
                "--ty", f"{ty:.2f}",
                "--p0", str(plate),
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            output = result.stdout

            match = re.search(r"zScale\s*=\s*([-+]?\d+(?:\.\d+)?)", output)

            if match is None:
                raise ValueError(f"Could not find zScale in output:\n{output}")

            zScale = float(match.group(1))

            print(zScale)

            # zScale = os.system(f'./run.exe 1 {event} --x0 {x0} --y0 {y0} --tx {tx:.2f} --ty {ty:.2f} --p0 {plate}')
            # print('a', zScale)
            ntuple.Fill(event, energy, tx, ty, plate, zScale)
            sampling += 1
        if sampling >= 3:
            break

outfile.cd()
ntuple.Write()
outfile.Write()
outfile.Close()