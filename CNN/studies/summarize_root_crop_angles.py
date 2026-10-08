#!/usr/bin/env python3
"""Summarize and plot CNN/Poisson-fit angles from ROOT crop angle CSV files."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ANGLE_FIELDS = ("cnn_theta_mrad", "fit_theta_mrad")
ANGLE_CUTS_MRAD = (5.0, 6.0, 10.0, 15.0, 20.0)


def numeric(value: str) -> float:
    return float(value) if value.strip() else float("nan")


def read_rows(path: Path, score_min: float) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["cnn_presence_score"] = numeric(row["cnn_presence_score"])
        for field in ANGLE_FIELDS:
            row[field] = numeric(row[field])
        row["fit_quality"] = row["fit_quality"].strip()
    return [row for row in rows if float(row["cnn_presence_score"]) >= score_min]


def finite(rows: list[dict[str, object]], field: str) -> np.ndarray:
    return np.asarray([
        float(row[field]) for row in rows if np.isfinite(float(row[field]))
    ])


def quantiles(values: np.ndarray) -> dict[str, float]:
    if not values.size:
        return {}
    return {str(q): float(np.quantile(values, q)) for q in (0.0, 0.5, 0.9, 0.95, 0.99, 1.0)}


def summarize(name: str, rows: list[dict[str, object]]) -> dict[str, object]:
    coherent = [row for row in rows if row["fit_quality"] == "coherent"]
    cnn = finite(rows, "cnn_theta_mrad")
    fit = finite(coherent, "fit_theta_mrad")
    return {
        "name": name,
        "crops_after_score_filter": len(rows),
        "fit_coherent": len(coherent),
        "cnn_theta_mrad_quantiles": quantiles(cnn),
        "fit_theta_mrad_quantiles_coherent": quantiles(fit),
        "cnn_theta_mrad_at_least": {
            str(cut): int(np.count_nonzero(cnn >= cut)) for cut in ANGLE_CUTS_MRAD
        },
        "fit_theta_mrad_at_least_coherent": {
            str(cut): int(np.count_nonzero(fit >= cut)) for cut in ANGLE_CUTS_MRAD
        },
    }


def plot_distributions(datasets: list[tuple[str, list[dict[str, object]]]], output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    bins = np.linspace(0.0, 40.0, 81)
    for name, rows in datasets:
        cnn = finite(rows, "cnn_theta_mrad")
        coherent = [row for row in rows if row["fit_quality"] == "coherent"]
        fit = finite(coherent, "fit_theta_mrad")
        axes[0].hist(cnn, bins=bins, histtype="step", linewidth=1.8, density=True, label=name)
        axes[1].hist(fit, bins=bins, histtype="step", linewidth=1.8, density=True, label=name)
    for axis, title in zip(axes, ("CNN angle", "Poisson-hot fit angle (coherent only)"), strict=True):
        axis.set(xlabel="theta [mrad]", ylabel="density", title=title, xlim=(0, 40))
        axis.grid(alpha=0.2)
        axis.legend()
    figure.savefig(output, dpi=160)
    plt.close(figure)


def plot_cnn_vs_fit(datasets: list[tuple[str, list[dict[str, object]]]], output: Path) -> None:
    figure, axis = plt.subplots(figsize=(6.2, 5.4), constrained_layout=True)
    for name, rows in datasets:
        coherent = [row for row in rows if row["fit_quality"] == "coherent"]
        cnn = finite(coherent, "cnn_theta_mrad")
        fit = finite(coherent, "fit_theta_mrad")
        axis.scatter(cnn, fit, s=7, alpha=0.3, label=name)
    axis.plot([0, 40], [0, 40], color="black", linestyle="--", linewidth=1, label="CNN = fit")
    axis.set(
        xlabel="theta CNN [mrad]", ylabel="theta fit [mrad]",
        title="CNN vs Poisson-hot fit", xlim=(0, 40), ylim=(0, 40), aspect="equal",
    )
    axis.grid(alpha=0.2)
    axis.legend()
    figure.savefig(output, dpi=160)
    plt.close(figure)


def write_high_angle_rows(datasets: list[tuple[str, list[dict[str, object]]]], output: Path, count: int) -> None:
    rows = []
    for name, source_rows in datasets:
        for row in source_rows:
            rows.append({
                "source": name,
                "root_entry": row["root_entry"],
                "cell_id": row["cell_id"],
                "tag_cell_id": row["tag_cell_id"],
                "cnn_presence_score": row["cnn_presence_score"],
                "cnn_theta_mrad": row["cnn_theta_mrad"],
                "fit_theta_mrad": row["fit_theta_mrad"],
                "fit_quality": row["fit_quality"],
                "cnn_window_x_start": row["cnn_window_x_start"],
                "cnn_window_y_start": row["cnn_window_y_start"],
            })
    rows.sort(key=lambda row: float(row["cnn_theta_mrad"]), reverse=True)
    fields = list(rows[0]) if rows else ["source", "root_entry"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows[:count])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--score-min", type=float, default=0.0)
    parser.add_argument("--top", type=int, default=100)
    args = parser.parse_args()
    if not 0.0 <= args.score_min <= 1.0:
        parser.error("--score-min must be in [0, 1]")
    if args.top <= 0:
        parser.error("--top must be positive")

    datasets = [(path.stem.removesuffix("_angles"), read_rows(path, args.score_min)) for path in args.input]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_distributions(datasets, args.output_dir / "angle_distributions.png")
    plot_cnn_vs_fit(datasets, args.output_dir / "cnn_vs_fit.png")
    write_high_angle_rows(datasets, args.output_dir / "top_cnn_angles.csv", args.top)
    summary = {
        "score_min": args.score_min,
        "angle_histogram_range_mrad": [0.0, 40.0],
        "files": [summarize(name, rows) for name, rows in datasets],
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir.resolve()), **summary}, indent=2))


if __name__ == "__main__":
    main()
