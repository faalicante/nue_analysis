import ROOT
import re
import os
from itertools import product
from collections import Counter

ROOT.gROOT.SetBatch(True)

def InterpolateBadBins(h, dipFrac = 0.2, neighCompat = 0.5):
    nb = h.GetNbinsX()
    for i in range(1, nb):
        yL = h.GetBinContent(i-1)
        y0 = h.GetBinContent(i)
        yR = h.GetBinContent(i+1)
        if (yR <= dipFrac * yL): yR = h.GetBinContent(i+2)
        if (yR <= dipFrac * yL): yR = h.GetBinContent(i+3)
        if (yL <= 0 or yR <= dipFrac * yL): continue
        
        isLow = (y0 < dipFrac * yL) and (y0 < dipFrac * yR)
        eps = 1e-9
        compatibleNeighbors = abs(yL - yR) / (max(yL, yR) + eps) < neighCompat
        # print(isLow, compatibleNeighbors, i, yL, y0, yR, abs(yL - yR) / (max(yL, yR) + eps))
        if (isLow and compatibleNeighbors):
            yInterp = 0.5 * (yL + yR)
            h.SetBinContent(i, yInterp)

path = '/eos/experiment/sndlhc/users/falicant'
opath = '/eos/user/f/falicant/shift/long_spectra_2'

file = open('doubles.txt', 'r')
lines = file.readlines()
file.close()


survey_sigma = [2,3]
survey_thr = [0.05]
mods = ['original', 'interpolated', 'smoothed', 'smoothed50']

current_folder = None
results = dict()

for line in lines:
    line = line.strip()
    if line.startswith('RUN'):
        current_folder = line
    elif line and current_folder:
        c = ROOT.TCanvas('c', 'c', 1500, 1800)
        c.Divide(2, 4)
        parts = re.split(r'\s*\*\s*', line)
        if len(parts) >= 3:
            # if current_folder != 'RUN1/b121' or parts[1] != '22': continue
            file_path = current_folder + f'/shift/peaks/peaks_{parts[1]}_p.root'
            if not os.path.exists(file_path):
                print(f"File not found: {file_path}")
                continue 
            print(f"Processing file: {file_path}")
            peak_file = ROOT.TFile.Open(file_path)
            h = peak_file.Get(f'h_long_{parts[2]}_{parts[0]}')
            h.SetDirectory(0)
            # InterpolateBadBins(h)
            # h.Smooth()
            s = ROOT.TSpectrum(5)
            hs = []
            for i, (mod, sigma, thr) in enumerate(product(mods,survey_sigma, survey_thr)):
                hs.append(h.Clone(f'h_{sigma}_{thr}'))
                hs[i].SetDirectory(0)
                if mod == 'interpolated':
                    InterpolateBadBins(hs[i])
                if mod == 'smoothed': 
                    InterpolateBadBins(hs[i])
                    hs[i].Smooth()
                if mod == 'smoothed50':
                    InterpolateBadBins(hs[i])
                    hs[i].Smooth(50)
                c.cd(i+1)
                hs[i].SetTitle(f'sigma={sigma}, thr={thr}, {mod}') 
                hs[i].Draw()
                search_result = s.Search(hs[i], sigma, 'nobackground', thr)
                results.setdefault((sigma,thr,mod), []).append(search_result)
            c.Update()
            c.SaveAs(opath+f'/{current_folder.split("/")[0]}_{current_folder.split("/")[1]}_{parts[1]}_{parts[2]}.png')        
            peak_file.Close()
            del c

# Now count frequencies for each (sigma, thr) pair
for key, values in results.items():
    if isinstance(key, tuple):  # Only process (sigma, thr) keys
        counts = Counter(values)
        print(f"For sigma={key[0]}, thr={key[1]}, mod={key[2]}:")
        for result, count in counts.items():
            print(f"  Result {result}: {count} times")