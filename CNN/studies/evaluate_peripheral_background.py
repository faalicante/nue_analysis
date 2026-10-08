#!/usr/bin/env python3
"""Calibrate scan scores using signal-centred maps away from the MC trajectory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scanning.root_xypseg_to_scan_hdf5 import coordinate_to_local_bin


def distance_from_window_to_trajectory(
    x_start: int,
    y_start: int,
    crop_size: int,
    trajectory_x: np.ndarray,
    trajectory_y: np.ndarray,
) -> float:
    """Minimum Euclidean bin distance between a crop rectangle and an MC path."""
    x_stop = x_start + crop_size
    y_stop = y_start + crop_size
    dx = np.maximum.reduce((x_start - trajectory_x, np.zeros_like(trajectory_x), trajectory_x - x_stop))
    dy = np.maximum.reduce((y_start - trajectory_y, np.zeros_like(trajectory_y), trajectory_y - y_stop))
    return float(np.min(np.hypot(dx, dy)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--volumes", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--far-margin-bins", type=float, default=15.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    peripheral_scores: list[float] = []
    signal_scores: list[float] = []
    ignored_scores: list[float] = []
    event_signal_maxima: list[float] = []
    event_peripheral_maxima: list[float] = []
    event_rows: list[dict[str, object]] = []

    with h5py.File(args.volumes, "r") as source, h5py.File(args.predictions, "r") as prediction:
        crop_size = int(prediction.attrs["crop_size"])
        x_positions = prediction["grid_x_start_bin"][:]
        y_positions = prediction["grid_y_start_bin"][:]
        for index, event_id in enumerate(source["event_id"][:]):
            score_grid = prediction["window_presence_score"][index]
            x_edges = source["x_edges_um"][index]
            y_edges = source["y_edges_um"][index]
            vertex_x = coordinate_to_local_bin(x_edges, float(source["truth/xpos_um"][index]))
            vertex_y = coordinate_to_local_bin(y_edges, float(source["truth/ypos_um"][index]))
            slope_x, slope_y = source["truth/slope_xy_bins_per_z"][index]
            crop_plate = int(source["truth/plate"][index]) - 3
            development = min(40, 57 - crop_plate)
            z_offset = np.arange(development + 1, dtype=np.float64)
            trajectory_x = vertex_x + float(slope_x) * z_offset
            trajectory_y = vertex_y + float(slope_y) * z_offset
            target_x = float(source["truth/propagated_center_x_bin"][index])
            target_y = float(source["truth/propagated_center_y_bin"][index])

            counts = {"signal": 0, "peripheral": 0, "ignored": 0}
            event_signal_scores: list[float] = []
            event_peripheral_scores: list[float] = []
            for row, y_start in enumerate(y_positions):
                for col, x_start in enumerate(x_positions):
                    score = float(score_grid[row, col])
                    contains_target = (
                        x_start <= target_x < x_start + crop_size
                        and y_start <= target_y < y_start + crop_size
                    )
                    distance = distance_from_window_to_trajectory(
                        int(x_start), int(y_start), crop_size, trajectory_x, trajectory_y
                    )
                    if contains_target:
                        signal_scores.append(score)
                        event_signal_scores.append(score)
                        counts["signal"] += 1
                    elif distance >= args.far_margin_bins:
                        peripheral_scores.append(score)
                        event_peripheral_scores.append(score)
                        counts["peripheral"] += 1
                    else:
                        ignored_scores.append(score)
                        counts["ignored"] += 1
            signal_maximum = max(event_signal_scores)
            peripheral_maximum = max(event_peripheral_scores)
            event_signal_maxima.append(signal_maximum)
            event_peripheral_maxima.append(peripheral_maximum)
            event_rows.append({
                "event_id": int(event_id), **counts,
                "maximum_signal_window_score": signal_maximum,
                "maximum_peripheral_window_score": peripheral_maximum,
            })

    peripheral = np.asarray(peripheral_scores, dtype=np.float64)
    signal = np.asarray(signal_scores, dtype=np.float64)
    event_signal_max = np.asarray(event_signal_maxima, dtype=np.float64)
    event_peripheral_max = np.asarray(event_peripheral_maxima, dtype=np.float64)
    if peripheral.size == 0 or signal.size == 0:
        raise RuntimeError("empty signal or peripheral-background selection")

    thresholds = [0.208049014210701, 0.5, 0.9, 0.95, 0.98, 0.99, 0.995]
    operating_points = [
        {
            "threshold": threshold,
            "signal_window_efficiency": float(np.mean(signal >= threshold)),
            "peripheral_window_fpr": float(np.mean(peripheral >= threshold)),
            "peripheral_windows_above": int(np.sum(peripheral >= threshold)),
            "signal_event_efficiency": float(np.mean(event_signal_max >= threshold)),
            "events_with_peripheral_candidate_fraction": float(
                np.mean(event_peripheral_max >= threshold)
            ),
        }
        for threshold in thresholds
    ]
    zero_observed_fp_threshold = float(np.nextafter(peripheral.max(), 1.0))
    report = {
        "volumes": str(args.volumes),
        "predictions": str(args.predictions),
        "far_margin_bins": args.far_margin_bins,
        "far_margin_um_nominal": args.far_margin_bins * 50.0,
        "events": event_rows,
        "counts": {
            "signal_windows": int(signal.size),
            "peripheral_background_windows": int(peripheral.size),
            "ignored_near_trajectory_windows": int(len(ignored_scores)),
        },
        "score_quantiles": {
            "signal": {str(q): float(np.quantile(signal, q)) for q in (0, 0.1, 0.5, 0.9, 1)},
            "peripheral": {str(q): float(np.quantile(peripheral, q)) for q in (0, 0.5, 0.9, 0.95, 0.99, 1)},
        },
        "operating_points": operating_points,
        "zero_observed_fp_threshold": zero_observed_fp_threshold,
        "signal_window_efficiency_at_zero_observed_fp": float(np.mean(signal >= zero_observed_fp_threshold)),
        "warning": "Correlated windows from three signal ROOTs; not a replacement for an independent hard-background sample.",
    }
    report_path = args.output_dir / "peripheral_background_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    bins = np.linspace(0.0, 1.0, 51)
    axis.hist(peripheral, bins=bins, density=True, histtype="step", linewidth=2,
              label=f"fondo periferico (n={peripheral.size})")
    axis.hist(signal, bins=bins, density=True, histtype="step", linewidth=2,
              label=f"finestre sul segnale (n={signal.size})")
    axis.set(xlabel="presence score", ylabel="densità", yscale="log", ylim=(1e-2, None))
    axis.grid(alpha=0.25)
    axis.legend()
    plot_path = args.output_dir / "peripheral_vs_signal_scores.png"
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    print(json.dumps({**report, "plot": str(plot_path), "report": str(report_path)}, indent=2))


if __name__ == "__main__":
    main()
