brickIDs=(12 13 14 21 22 23 24 31 32 33 34 41 42 43 44 51 52 53 54)
for ibrick in $(seq 0 18)
    do 
        cd b0000${brickIDs[ibrick]}
        echo "Processing nue brick ${brickIDs[ibrick]}"
        root -l hist_XYP_nue.root /afs/cern.ch/work/f/falicant/public/nue_analysis/tag_basetrack.C
        cd ../
    done

for ibrick in 0 4 8 12 16
    do
        cd b0000${brickIDs[ibrick]}
        echo "Processing muon brick ${brickIDs[ibrick]}"
        root -l hist_XYP_muon.root /afs/cern.ch/work/f/falicant/public/nue_analysis/tag_basetrack.C
        cd ../
    done
echo finish