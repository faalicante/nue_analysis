#!/usr/bin/env python3
"""Plot Poisson-hot XZ/YZ projections for selected scan-candidate ranks."""

from __future__ import annotations

import argparse
import base64
import json
import math
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import PowerNorm

from visualization.make_background_candidate_gifs import poisson_hot_threshold


Z_STEP_UM = 1350.0


def line_from_slope(start_um: float, slope_mrad: float, start_layer: int) -> np.ndarray:
    layers = np.arange(1, 58, dtype=float)
    return start_um + slope_mrad * Z_STEP_UM * (layers - start_layer) / 1000.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ranks", type=int, nargs="+", required=True)
    parser.add_argument("--output-png", type=Path, required=True)
    parser.add_argument("--output-html", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    by_rank = {int(row["rank"]): row for row in manifest["candidates"]}
    rows = [by_rank[rank] for rank in args.ranks]
    figure, axes = plt.subplots(
        len(rows), 2, figsize=(12.5, 4.4 * len(rows)), constrained_layout=True
    )
    axes = np.atleast_2d(axes)

    for row_index, candidate in enumerate(rows):
        with h5py.File(candidate["volume_path"], "r") as source:
            event_index = int(candidate["event_index"])
            x0 = int(candidate["display_x_start"])
            y0 = int(candidate["display_y_start"])
            size = int(candidate["display_crop_size"])
            raw = np.asarray(
                source["volumes_raw"][event_index, :, y0:y0 + size, x0:x0 + size],
                dtype=float,
            )
            mu = float(source["background_mu"][event_index])
            x_edges = np.asarray(source["x_edges_um"][event_index, x0:x0 + size + 1])
            y_edges = np.asarray(source["y_edges_um"][event_index, y0:y0 + size + 1])

        threshold = poisson_hot_threshold(mu)
        significant = np.where(raw >= threshold, np.maximum(raw - mu, 0.0), 0.0)
        projections = (significant.sum(axis=1), significant.sum(axis=2))
        edges = (x_edges / 1000.0, y_edges / 1000.0)
        coordinate_names = ("x", "y")
        cnn_slopes = (
            float(candidate["tx_pred_mrad"]), float(candidate["ty_pred_mrad"])
        )
        bin_widths = (float(np.median(np.diff(x_edges))), float(np.median(np.diff(y_edges))))
        component_slopes = (
            1000.0 * float(candidate["selected_component_slope_x"]) * bin_widths[0] / Z_STEP_UM,
            1000.0 * float(candidate["selected_component_slope_y"]) * bin_widths[1] / Z_STEP_UM,
        )
        starts = (float(candidate["start_x_um"]), float(candidate["start_y_um"]))
        start_layer = int(candidate["start_layer"])

        for column, (projection, physical_edges, coordinate, cnn_slope, component_slope, start) in enumerate(
            zip(projections, edges, coordinate_names, cnn_slopes, component_slopes, starts, strict=True)
        ):
            axis = axes[row_index, column]
            positive = projection[projection > 0]
            vmax = float(np.quantile(positive, 0.995)) if positive.size else 1.0
            mesh = axis.pcolormesh(
                physical_edges,
                np.arange(0.5, 58.5, 1.0),
                projection,
                cmap="viridis",
                norm=PowerNorm(gamma=0.55, vmin=0.0, vmax=max(vmax, 1.0)),
                shading="flat",
            )
            layers = np.arange(1, 58)
            cnn_line = line_from_slope(start, cnn_slope, start_layer) / 1000.0
            component_line = line_from_slope(start, component_slope, start_layer) / 1000.0
            axis.plot(cnn_line, layers, color="#ff3b30", lw=2.0, label=f"CNN: {cnn_slope:+.1f} mrad")
            if np.isfinite(component_slope):
                axis.plot(
                    component_line, layers, color="#ffffff", lw=1.8, ls="--",
                    label=f"fit hot: {component_slope:+.1f} mrad",
                )
            axis.scatter([start / 1000.0], [start_layer], s=42, marker="x", color="#ffcc00", lw=2.2, label="inizio stimato")
            axis.set_xlim(physical_edges[0], physical_edges[-1])
            axis.set_ylim(0.5, 57.5)
            axis.set_xlabel(f"{coordinate} [mm]")
            axis.set_ylabel("layer XYPseg")
            axis.set_title(f"Rank {candidate['rank']} — {coordinate.upper()}Z")
            axis.grid(alpha=0.16)
            axis.legend(loc="best", fontsize=9, framealpha=0.82)
            figure.colorbar(mesh, ax=axis, pad=0.01, label="eccesso q-hot integrato")

        cnn_theta = float(candidate["theta_pred_mrad"])
        component_theta = 1000.0 * math.atan(math.hypot(*component_slopes) / 1000.0)
        axes[row_index, 0].text(
            0.01, 1.13,
            f"Rank {candidate['rank']} · score={float(candidate['presence_score']):.3f} · "
            f"θ CNN={cnn_theta:.1f} mrad · θ fit hot={component_theta:.1f} mrad",
            transform=axes[row_index, 0].transAxes, fontsize=11, weight="bold",
        )

    figure.suptitle(
        "Proiezioni longitudinali dei candidati — voxel sopra soglia Poisson",
        fontsize=15,
    )
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_png, dpi=150, bbox_inches="tight")
    plt.close(figure)

    encoded = base64.b64encode(args.output_png.read_bytes()).decode("ascii")
    fragment = f'''<div id="candidate-xzyz-ranks-1-2" style="width:100%; color:var(--foreground);">
  <img src="data:image/png;base64,{encoded}" alt="Proiezioni XZ e YZ dei candidati rank 1 e rank 2 con slope CNN e fit geometrico" style="display:block; width:100%; height:auto;" />
</div>
'''
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    args.output_html.write_text(fragment, encoding="utf-8")
    print(json.dumps({
        "png": str(args.output_png.resolve()),
        "html": str(args.output_html.resolve()),
        "ranks": args.ranks,
    }, indent=2))


if __name__ == "__main__":
    main()
