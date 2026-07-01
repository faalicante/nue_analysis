import ROOT

tag = 'tagged_showers.txt'
sel = 'candidate_showers.txt'
stepZ      = 1350
tol = 1000

print("Cell, Combination, Tag, X, Y, Tx, Ty, Vx, Vy, Vidx")
with open(tag, 'r') as f1:
    lines = f1.readlines()
    for iline, line in enumerate(lines):
        if iline < 3: continue # Skip header lines
        line = line.replace('\n', '')
        line = line.replace('*', '')
        line = line.split()
        cell = int(line[1])
        comb = int(line[2])
        tag = int(line[3])
        x = int(line[4])
        y = int(line[5])
        tx = int(line[6])
        ty = int(line[7])
        with open(sel, 'r') as f2:
            sel_lines = f2.readlines()
            for jline, sel_line in enumerate(sel_lines):
                if jline < 3: continue
                sel_line = sel_line.replace('\n', '')
                sel_line = sel_line.replace('*', '')
                sel_line = sel_line.split()
                vcell = int(sel_line[2])
                vx = float(sel_line[3])
                vy = float(sel_line[4])
                plate = int(sel_line[6])
                vidx = int(sel_line[5])
                shiftX = tx / 1000.0 * stepZ * (plate - 1)
                shiftY = ty / 1000.0 * stepZ * (plate - 1)
                if vcell == cell and abs(vx - (x - shiftX)) < tol and abs(vy - (y - shiftY)) < tol:
                    print(f'cell: {cell}, comb: {comb}, tx: {tx}, ty: {ty}, tag: {tag}, x: {x}, y: {y}, vx: {int(vx)}, vy: {int(vy)}, plate: {plate}, vidx: {vidx}')