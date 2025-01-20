brickIDs=(75 76 77 78 58 59 60 147 115 133 116 117 118 119 120 121 122 134 135 136 140 82 83 84 85 86 100 101 102 103 104 153 154 158)
for ibrick in $(seq 1 34)
 do 
	#echo "cell ${brickIDs[ibrick]}"
	python shift_map_data.py --cell ${brickIDs[ibrick]}
 done
