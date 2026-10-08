from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


OUTPUT = Path("/Users/fabioali/SND@LHC/nue_analysis/CNN/output/hmu_review_en")
summary_path = OUTPUT / "bkg_mu_by_brick_summary.csv"

with summary_path.open(newline="") as handle:
    rows = list(csv.DictReader(handle))
rows.sort(key=lambda row: int(row["brick"].removeprefix("b")))

colors = {"fitlt5_mu_residual": "#126782", "gap5_10_mu_high": "#D16C34"}
labels = {"fitlt5_mu_residual": "H − μ run", "gap5_10_mu_high": "gap5_10_mu_high run"}
x = np.arange(len(rows))
means = [float(row["mean_bkg_mu"]) for row in rows]

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.facecolor": "white",
    }
)
fig, axis = plt.subplots(figsize=(16, 8.6), facecolor="white")
fig.subplots_adjust(left=0.075, right=0.98, bottom=0.22, top=0.76)
ink = "#102A43"
muted = "#526777"
fig.suptitle("Mean background μ by brick", x=0.055, y=0.95, ha="left", fontsize=25, weight="bold", color=ink)
fig.text(
    0.055,
    0.885,
    "All valid data-scan cells in bricks with converted candidate TXT files",
    fontsize=16,
    color="#126782",
)
axis.bar(x, means, color=[colors[row["model"]] for row in rows], width=0.72)
axis.set_xticks(x, [row["brick"] for row in rows], rotation=55, ha="right")
axis.set_ylabel("Mean bkg_mu per cell", fontsize=15)
axis.set_xlabel("Brick", fontsize=15)
axis.set_ylim(0, max(means) * 1.16)
axis.grid(axis="y", alpha=0.22)
axis.legend(
    handles=[Patch(facecolor=colors[name], label=labels[name]) for name in colors],
    frameon=False,
    loc="upper left",
)
for xi, mean in zip(x, means):
    axis.text(xi, mean + max(means) * 0.017, f"{mean:.1f}", ha="center", va="bottom", fontsize=8.5, color=ink)
fig.text(
    0.055,
    0.055,
    "Colour identifies the model run associated with the converted candidate TXT file. The mean uses the scalar background_mu stored for each cell.",
    fontsize=12.5,
    color=muted,
)
fig.text(
    0.055,
    0.020,
    "The two non-HDF5 conflicted copies found in b000431 and b000541 were excluded; all other listed HDF5 files were included.",
    fontsize=12.5,
    color=muted,
)
for extension in ("png", "svg", "pdf"):
    fig.savefig(OUTPUT / f"mean_bkg_mu_by_brick_EN.{extension}", dpi=180)
plt.close(fig)
print("Wrote mean_bkg_mu_by_brick_EN.{png,svg,pdf}")
