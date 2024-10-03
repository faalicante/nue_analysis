brickIDs=(11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44 51 52 53 54)
for ibrick in $(seq 0 19)
    do 
        cd b0000${brickIDs[ibrick]}
        echo "Processing nue brick ${brickIDs[ibrick]}"
        python /afs/cern.ch/work/f/falicant/public/nue_analysis/genXYP_mc.py -b ${brickIDs[ibrick]}
        cd ../
    done

# for ibrick in 0 4 8 12 16
#     do
#         echo "Processing muon brick ${brickIDs[ibrick]}"
#         python /afs/cern.ch/work/f/falicant/public/nue_analysis/spherocity.py -b ${brickIDs[ibrick]}
#     done
# echo finish