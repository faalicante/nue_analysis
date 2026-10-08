from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


WORKSPACE = Path("/Users/fabioali/SND@LHC/nue_analysis/CNN")
OUTPUT = WORKSPACE / "output" / "hmu_review_en"
VOLUME_ROOT = Path("/Users/fabioali/cernbox/CNN/data_scan_volumes")
COUNTS = OUTPUT / "candidate_counts_by_brick.csv"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.facecolor": "white",
    }
)


def load_brick_mu(brick: str) -> tuple[np.ndarray, int, list[str]]:
    """Read scalar background_mu values once per stored data-scan cell."""
    values: list[np.ndarray] = []
    valid_files = 0
    skipped: list[str] = []
    for path in sorted((VOLUME_ROOT / brick).glob("*.h5")):
        if not h5py.is_hdf5(path):
            skipped.append(path.name)
            continue
        with h5py.File(path, "r") as source:
            if "background_mu" not in source:
                skipped.append(path.name)
                continue
            values.append(np.asarray(source["background_mu"], dtype=float).ravel())
            valid_files += 1
    if not values:
        raise RuntimeError(f"No valid background_mu values found for {brick}")
    return np.concatenate(values), valid_files, skipped


records: list[dict[str, object]] = []
grouped: dict[str, list[tuple[str, np.ndarray]]] = defaultdict(list)
with COUNTS.open(newline="") as handle:
    for row in csv.DictReader(handle):
        brick = row["brick"]
        values, valid_files, skipped = load_brick_mu(brick)
        model = row["model"]
        grouped[model].append((brick, values))
        records.append(
            {
                "brick": brick,
                "model": model,
                "valid_hdf5_files": valid_files,
                "cells": len(values),
                "min_bkg_mu": float(np.min(values)),
                "mean_bkg_mu": float(np.mean(values)),
                "p05_bkg_mu": float(np.quantile(values, 0.05)),
                "median_bkg_mu": float(np.median(values)),
                "p95_bkg_mu": float(np.quantile(values, 0.95)),
                "max_bkg_mu": float(np.max(values)),
                "skipped_non_hdf5_files": "; ".join(skipped),
            }
        )

summary_path = OUTPUT / "bkg_mu_by_brick_summary.csv"
with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)

all_values = np.concatenate([values for _, values in sum(grouped.values(), [])])
step = 0.25
left = max(0.0, np.floor(all_values.min() / step) * step)
right = np.ceil(all_values.max() / step) * step
edges = np.arange(left, right + step, step)

ink = "#102A43"
blue = "#126782"
muted = "#526777"
fig, axes = plt.subplots(1, 2, figsize=(16, 9.2), sharey=True, facecolor="white")
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.38, top=0.76, wspace=0.13)
fig.suptitle("Background μ distributions by brick", x=0.055, y=0.955, ha="left", fontsize=25, weight="bold", color=ink)
fig.text(
    0.055,
    0.895,
    "All valid data-scan cells in bricks with converted candidate TXT files",
    fontsize=16,
    color=blue,
)

labels = {
    "fitlt5_mu_residual": "H − μ candidate TXT files (18 bricks)",
    "gap5_10_mu_high": "gap5_10_mu_high candidate TXT files (6 bricks)",
}
for axis, model in zip(axes, ["fitlt5_mu_residual", "gap5_10_mu_high"]):
    items = grouped[model]
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(items)))
    for (brick, values), color in zip(items, colors):
        axis.hist(
            values,
            bins=edges,
            density=True,
            histtype="step",
            linewidth=1.65,
            color=color,
            label=f"{brick}  (n={len(values)})",
        )
    axis.set_title(labels[model], loc="left", fontsize=14, color=ink, pad=16)
    axis.set_xlabel("bkg_mu per cell", fontsize=15)
    axis.grid(axis="y", alpha=0.20)
    axis.set_xlim(left, right)
    columns = 3 if len(items) > 9 else 2
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.20),
        ncol=columns,
        frameon=False,
        fontsize=8.5,
        columnspacing=1.35,
        handlelength=1.7,
    )

axes[0].set_ylabel("Normalized cell density", fontsize=15)
fig.text(
    0.055,
    0.095,
    "Each curve is normalized within its brick. Values come from the scalar background_mu field stored for each data-scan cell.",
    fontsize=12.5,
    color=muted,
)
fig.text(
    0.055,
    0.052,
    "The two non-HDF5 conflicted copies found in b000431 and b000541 were excluded; all other listed HDF5 files were included.",
    fontsize=12.5,
    color=muted,
)

for extension in ("png", "svg", "pdf"):
    fig.savefig(OUTPUT / f"bkg_mu_by_brick_overlay_EN.{extension}", dpi=180)
plt.close(fig)
print(f"Wrote {summary_path}")
