#!/usr/bin/env python3
"""Render XZ/YZ projections of ROOT crops selected from an angle CSV."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
import numpy as np
import uproot
from matplotlib.colors import PowerNorm

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from visualization.make_background_candidate_gifs import poisson_hot_threshold, root_th2_smooth


def number(row: dict[str, str], field: str) -> float:
    value = row.get(field, "").strip()
    return float(value) if value else float("nan")


def selected_rows(
    path: Path,
    *,
    angle_field: str,
    min_angle: float,
    score_min: float,
    coherent_fit_only: bool,
    limit: int | None,
) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = [
        row for row in rows
        if math.isfinite(number(row, angle_field))
        and number(row, angle_field) >= min_angle
        and number(row, "cnn_presence_score") >= score_min
        and (not coherent_fit_only or row.get("fit_quality") == "coherent")
    ]
    selected.sort(key=lambda row: number(row, angle_field), reverse=True)
    return selected[:limit] if limit is not None else selected


def slope_line(start: float, slope: float, start_layer: int) -> np.ndarray:
    return start + slope * (np.arange(1, 58) - start_layer)


def render(
    raw: np.ndarray,
    background_mu: float,
    row: dict[str, str],
    output: Path,
    *,
    smooth: bool,
    poisson_tail_probability: float,
) -> dict[str, object]:
    values = root_th2_smooth(raw) if smooth else raw.astype(float, copy=False)
    threshold = poisson_hot_threshold(background_mu, poisson_tail_probability)
    significant = np.where(values >= threshold, np.maximum(values - background_mu, 0.0), 0.0)
    projections = (significant.sum(axis=1), significant.sum(axis=2))
    cnn_slopes = (number(row, "cnn_slope_x_bins_per_layer"), number(row, "cnn_slope_y_bins_per_layer"))
    fit_slopes = (number(row, "fit_slope_x_bins_per_layer"), number(row, "fit_slope_y_bins_per_layer"))
    start_layer = number(row, "fit_start_layer")
    start_layer = int(start_layer) if math.isfinite(start_layer) else 1
    starts = (number(row, "fit_start_x_bin"), number(row, "fit_start_y_bin"))
    starts = tuple(value if math.isfinite(value) else 16.0 for value in starts)

    figure, axes = plt.subplots(1, 2, figsize=(12.4, 5.2), constrained_layout=True)
    for axis, projection, coordinate, cnn_slope, fit_slope, start in zip(
        axes, projections, ("x", "y"), cnn_slopes, fit_slopes, starts, strict=True
    ):
        positive = projection[projection > 0]
        vmax = float(np.quantile(positive, 0.995)) if positive.size else 1.0
        mesh = axis.pcolormesh(
            np.arange(0.5, 58.5), np.arange(0.0, 33.0), projection.T,
            cmap="viridis", norm=PowerNorm(gamma=0.55, vmin=0.0, vmax=max(vmax, 1.0)),
            shading="flat",
        )
        layers = np.arange(1, 58)
        axis.plot(layers, slope_line(start, cnn_slope, start_layer), color="#ff3b30", lw=2,
                  label=f"CNN {cnn_slope:+.2f} bin/layer")
        if all(math.isfinite(value) for value in fit_slopes):
            axis.plot(layers, slope_line(start, fit_slope, start_layer), color="white", lw=1.8,
                      ls="--", label=f"fit {fit_slope:+.2f} bin/layer")
        axis.scatter([start_layer], [start], color="#ffcc00", marker="x", s=44, lw=2.2,
                     label="inizio fit")
        axis.set(xlim=(0.5, 57.5), ylim=(0.0, 32.0), xlabel="layer XYPseg (z)",
                 ylabel=f"{coordinate} [bin]", title=f"{coordinate.upper()}Z")
        axis.grid(alpha=0.16)
        axis.legend(loc="best", fontsize=8.5, framealpha=0.84)
        figure.colorbar(mesh, ax=axis, pad=0.01, label="eccesso q-hot integrato")
    figure.suptitle(
        f"ROOT entry {row['root_entry']} · cell {row['cell_id']} · tag {row['tag_cell_id']} · "
        f"score={number(row, 'cnn_presence_score'):.3f} · "
        f"theta CNN={number(row, 'cnn_theta_mrad'):.1f} mrad · "
        f"theta fit={number(row, 'fit_theta_mrad'):.1f} mrad · {row['fit_quality']}",
        fontsize=12.5,
    )
    figure.savefig(output, dpi=135, bbox_inches="tight")
    plt.close(figure)
    return {
        "root_entry": int(row["root_entry"]),
        "cell_id": int(row["cell_id"]),
        "tag_cell_id": int(row["tag_cell_id"]),
        "cnn_presence_score": number(row, "cnn_presence_score"),
        "cnn_theta_mrad": number(row, "cnn_theta_mrad"),
        "fit_theta_mrad": number(row, "fit_theta_mrad"),
        "fit_quality": row["fit_quality"],
        "poisson_hot_count_threshold": threshold,
        "projection_smoothed": smooth,
        "projection_poisson_tail_probability": poisson_tail_probability,
        "image": str(output.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--angles", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--angle-source", choices=("cnn", "fit"), default="cnn")
    parser.add_argument("--min-angle-mrad", type=float, default=8.0)
    parser.add_argument("--score-min", type=float, default=0.0)
    parser.add_argument("--coherent-fit-only", action="store_true")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--smooth", action="store_true")
    parser.add_argument("--poisson-tail-probability", type=float, default=1.0e-4)
    args = parser.parse_args()
    if args.min_angle_mrad < 0.0:
        parser.error("--min-angle-mrad must be non-negative")
    if not 0.0 <= args.score_min <= 1.0:
        parser.error("--score-min must be in [0, 1]")
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be positive")
    if not 0.0 < args.poisson_tail_probability < 1.0:
        parser.error("--poisson-tail-probability must be between 0 and 1")

    angle_field = f"{args.angle_source}_theta_mrad"
    rows = selected_rows(
        args.angles, angle_field=angle_field, min_angle=args.min_angle_mrad,
        score_min=args.score_min, coherent_fit_only=args.coherent_fit_only,
        limit=args.limit,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rendered: list[dict[str, object]] = []
    with uproot.open(args.input_root) as root_file:
        tree = root_file["samples"]
        for rank, row in enumerate(rows, start=1):
            entry = int(row["root_entry"])
            arrays = tree.arrays(
                ["counts", "background_mu"], entry_start=entry, entry_stop=entry + 1, library="np"
            )
            filename = f"rank_{rank:03d}_entry_{entry:05d}_theta_{number(row, angle_field):06.2f}_xz_yz.png"
            rendered.append({
                "rank": rank,
                "angle_source": args.angle_source,
                "selected_angle_mrad": number(row, angle_field),
                **render(
                    np.asarray(arrays["counts"][0], dtype=float), float(arrays["background_mu"][0]),
                    row, args.output_dir / filename, smooth=args.smooth,
                    poisson_tail_probability=args.poisson_tail_probability,
                ),
            })
    with (args.output_dir / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = list(rendered[0]) if rendered else ["rank", "root_entry"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rendered)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "input_root": str(args.input_root.resolve()), "angles": str(args.angles.resolve()),
        "angle_source": args.angle_source, "min_angle_mrad": args.min_angle_mrad,
        "score_min": args.score_min, "coherent_fit_only": args.coherent_fit_only,
        "count": len(rendered), "candidates": rendered,
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"count": len(rendered), "output_dir": str(args.output_dir.resolve())}, indent=2))


if __name__ == "__main__":
    main()
