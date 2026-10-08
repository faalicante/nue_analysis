#!/usr/bin/env python3
"""Compare the original and capped crop centres for signal event 938."""

from __future__ import annotations

from pathlib import Path

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import uproot


EVENT = 938
ROOT_PATH = Path("/Users/fabioali/cernbox/shift/nue_Euniform_FLUKA/b000021.0.0.939.trk.root")
OUTPUT = Path("runs/cnn21d_sample_type_no_poisson_theta10_50/pilot/truth_tracks/event_938_crop_before_after.png")
N_LAYERS = 57
DZ_UM = 1350.0
XY_BIN_UM = 50.0
SHRINK = 1315.0 / 1350.0


def main() -> None:
    # Metadata from nue_int_10k.txt for event 938.
    vertex_x, vertex_y = 199764.91, 186243.88
    tx_mrad, ty_mrad, interaction_layer = -49.57, -12.40, 0
    tx = tx_mrad * SHRINK / 1000.0
    ty = ty_mrad * SHRINK / 1000.0

    branches = ["s/s.eFlag", "s/s.eMCTrack", "s/s.eX", "s/s.eY", "s/s.eZ"]
    arrays = uproot.open(ROOT_PATH)["tracks"].arrays(branches, library="ak")
    selected = arrays[branches[0]] == 1
    flag, mc_track, x_um, y_um, z_um = (
        ak.to_numpy(ak.flatten(arrays[name][selected])) for name in branches
    )
    del flag
    layer = (z_um + DZ_UM * (N_LAYERS - 1)) / DZ_UM
    x_vertex_bins = (x_um - vertex_x) / XY_BIN_UM
    y_vertex_bins = (y_um - vertex_y) / XY_BIN_UM
    primary = mc_track == 1

    variants = [
        ("Prima: propagazione su 57 piatti", N_LAYERS - interaction_layer),
        ("Dopo: min(40, 57 − pn)", min(40, N_LAYERS - interaction_layer)),
    ]
    figure, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True, sharex="col", sharey="col")
    colors = layer

    for row, (title, propagated) in enumerate(variants):
        center_offset = 0.5 * propagated
        center_x_bins = tx * DZ_UM * center_offset / XY_BIN_UM
        center_y_bins = ty * DZ_UM * center_offset / XY_BIN_UM

        for axis, coordinate, center, coordinate_name in (
            (axes[row, 0], x_vertex_bins, center_x_bins, "x"),
            (axes[row, 1], y_vertex_bins, center_y_bins, "y"),
        ):
            axis.scatter(layer, coordinate, c=colors, cmap="viridis", s=10, alpha=0.45, linewidths=0)
            axis.scatter(layer[primary], coordinate[primary], c="crimson", s=35, label="eMCTrack=1")
            axis.axhspan(center - 16, center + 16, color="gold", alpha=0.10, label="crop 32")
            axis.axhspan(center - 10, center + 10, color="limegreen", alpha=0.14, label="crop 20")
            axis.axhline(center, color="black", linestyle=":", lw=1.5, label="centro")
            axis.set_ylabel(f"{coordinate_name} − vertice [bin]" if row == 0 else f"{coordinate_name} − vertice [bin]")
            axis.grid(alpha=0.18)

        xy_axis = axes[row, 2]
        xy_axis.scatter(x_vertex_bins, y_vertex_bins, c=colors, cmap="viridis", s=10, alpha=0.45, linewidths=0)
        xy_axis.scatter(x_vertex_bins[primary], y_vertex_bins[primary], c="crimson", s=38)
        xy_axis.add_patch(Rectangle((center_x_bins - 16, center_y_bins - 16), 32, 32, fill=False, ec="goldenrod", lw=2))
        xy_axis.add_patch(Rectangle((center_x_bins - 10, center_y_bins - 10), 20, 20, fill=False, ec="limegreen", lw=2))
        xy_axis.scatter([center_x_bins], [center_y_bins], marker="x", c="black", s=70, linewidths=2)
        xy_axis.set_aspect("equal", adjustable="box")
        xy_axis.set_ylabel("y − vertice [bin]")

        axes[row, 0].set_title(f"{title}\ncentro dopo {center_offset:g} piatti — proiezione z–x")
        axes[row, 1].set_title(f"Proiezione z–y")
        axes[row, 2].set_title("Vista xy")

    for axis in axes[:, 0:2].flat:
        axis.set_xlabel("layer z")
    for axis in axes[:, 2]:
        axis.set_xlabel("x − vertice [bin]")
    axes[0, 0].legend(fontsize=8, loc="best")
    figure.suptitle(
        "Evento 938 — confronto reale del centro del crop\n"
        "riga superiore: vecchia formula; riga inferiore: limite di 40 piatti",
        fontsize=15,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=180)
    plt.close(figure)
    print(OUTPUT)


if __name__ == "__main__":
    main()
