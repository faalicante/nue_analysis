from __future__ import annotations

import csv
from pathlib import Path


OUTPUT = Path("/Users/fabioali/SND@LHC/nue_analysis/CNN/output/hmu_review_en")
counts_path = OUTPUT / "candidate_counts_by_brick.csv"
background_path = OUTPUT / "bkg_mu_by_brick_summary.csv"
result_path = OUTPUT / "candidate_rates_by_brick.csv"

with counts_path.open(newline="") as handle:
    counts = {row["brick"]: row for row in csv.DictReader(handle)}
with background_path.open(newline="") as handle:
    background = {row["brick"]: row for row in csv.DictReader(handle)}

if set(counts) != set(background):
    raise RuntimeError("Candidate and HDF5 brick lists do not match")

rows: list[dict[str, object]] = []
for brick, count in counts.items():
    source = background[brick]
    candidates = int(count["candidates"])
    cells = int(source["cells"])
    rows.append(
        {
            "brick": brick,
            "model": count["model"],
            "candidates": candidates,
            "valid_hdf5_files": int(source["valid_hdf5_files"]),
            "cells": cells,
            "candidates_per_cell": candidates / cells,
            "candidates_per_100_cells": 100 * candidates / cells,
        }
    )

with result_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote {result_path}")
