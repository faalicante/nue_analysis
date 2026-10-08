#!/usr/bin/env python3
"""Evaluate CNN score thresholds before and after the shower-size post-filter."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import re
from pathlib import Path

import h5py
import numpy as np

from studies.analyze_shower_size_features import size_features
from training.cnn21d.metrics import slopes_to_angles
from scanning.scan_cnn21d_volumes import cluster_contains_truth, connected_components


def candidate_size(
    source: h5py.File,
    event_index: int,
    row: int,
    col: int,
    x_positions: np.ndarray,
    y_positions: np.ndarray,
    cache: dict[tuple[int, int, int], int],
) -> int:
    key = (event_index, row, col)
    if key not in cache:
        x0 = int(x_positions[col])
        y0 = int(y_positions[row])
        raw = source["volumes_raw"][event_index, :, y0:y0 + 20, x0:x0 + 20]
        cache[key] = int(size_features(raw, float(source["background_mu"][event_index]))[
            "largest_component_voxels"
        ])
    return cache[key]


def evaluate_background(
    volume_glob: str,
    prediction_glob: str,
    thresholds: list[float],
    min_component_voxels: int,
) -> dict[float, dict[str, int]]:
    sources = {path.name: path for path in map(Path, glob.glob(volume_glob))}
    result = {
        threshold: {
            "background_clusters": 0,
            "background_cells": 0,
            "background_clusters_after_size": 0,
            "background_cells_after_size": 0,
        }
        for threshold in thresholds
    }
    for prediction_path in sorted(map(Path, glob.glob(prediction_glob))):
        source_name = re.sub(r"_t\d+\.h5$", ".h5", prediction_path.name)
        with h5py.File(sources[source_name], "r") as source, h5py.File(prediction_path, "r") as pred:
            scores = np.asarray(pred["window_presence_score"])
            x_positions = np.asarray(pred["grid_x_start_bin"], dtype=int)
            y_positions = np.asarray(pred["grid_y_start_bin"], dtype=int)
            cache: dict[tuple[int, int, int], int] = {}
            for threshold in thresholds:
                cells_before: set[int] = set()
                cells_after: set[int] = set()
                for event_index in range(len(scores)):
                    components = connected_components(scores[event_index] >= threshold)
                    if components:
                        cells_before.add(int(source["cell_id"][event_index]))
                    result[threshold]["background_clusters"] += len(components)
                    for component in components:
                        row, col = max(component, key=lambda rc: float(scores[event_index][rc]))
                        if candidate_size(
                            source, event_index, row, col, x_positions, y_positions, cache
                        ) >= min_component_voxels:
                            result[threshold]["background_clusters_after_size"] += 1
                            cells_after.add(int(source["cell_id"][event_index]))
                result[threshold]["background_cells"] += len(cells_before)
                result[threshold]["background_cells_after_size"] += len(cells_after)
    return result


def evaluate_signal(
    volumes: Path,
    predictions: Path,
    thresholds: list[float],
    theta_min: float,
    min_component_voxels: int,
    truth_propagation_limit: int | None,
) -> tuple[int, dict[float, dict[str, int]]]:
    result: dict[float, dict[str, int]] = {}
    with h5py.File(volumes, "r") as source, h5py.File(predictions, "r") as pred:
        scores = np.asarray(pred["window_presence_score"])
        theta = slopes_to_angles(np.asarray(source["truth/slope_xy_bins_per_z"]))[0]
        selected = np.flatnonzero(theta > theta_min)
        if truth_propagation_limit is None:
            truth_x = np.asarray(source["truth/propagated_center_x_bin"])
            truth_y = np.asarray(source["truth/propagated_center_y_bin"])
        else:
            crop_plate = np.asarray(source["truth/plate"], dtype=np.int64) - 3
            propagated = np.minimum(truth_propagation_limit, 57 - crop_plate)
            factor = 1350.0 * 0.5 * propagated / 1000.0
            center_x_um = (
                np.asarray(source["truth/xpos_um"])
                + np.asarray(source["truth/ltx_mrad_corrected"]) * factor
            )
            center_y_um = (
                np.asarray(source["truth/ypos_um"])
                + np.asarray(source["truth/lty_mrad_corrected"]) * factor
            )
            truth_x = np.asarray([
                np.interp(value, 0.5 * (edges[:-1] + edges[1:]), np.arange(len(edges) - 1))
                for value, edges in zip(center_x_um, source["x_edges_um"][:], strict=True)
            ])
            truth_y = np.asarray([
                np.interp(value, 0.5 * (edges[:-1] + edges[1:]), np.arange(len(edges) - 1))
                for value, edges in zip(center_y_um, source["y_edges_um"][:], strict=True)
            ])
        x_positions = np.asarray(pred["grid_x_start_bin"], dtype=int)
        y_positions = np.asarray(pred["grid_y_start_bin"], dtype=int)
        crop_size = int(pred.attrs["crop_size"])
        cache: dict[tuple[int, int, int], int] = {}
        for threshold in thresholds:
            any_candidate = 0
            truth_matched = 0
            any_after_size = 0
            truth_after_size = 0
            for event_index in selected:
                components = connected_components(scores[event_index] >= threshold)
                any_candidate += int(bool(components))
                event_truth = False
                event_any_size = False
                event_truth_size = False
                for component in components:
                    matched = cluster_contains_truth(
                        component, float(truth_x[event_index]), float(truth_y[event_index]),
                        x_positions.tolist(), y_positions.tolist(), crop_size,
                    )
                    event_truth |= matched
                    row, col = max(component, key=lambda rc: float(scores[event_index][rc]))
                    passes_size = candidate_size(
                        source, int(event_index), row, col, x_positions, y_positions, cache
                    ) >= min_component_voxels
                    event_any_size |= passes_size
                    event_truth_size |= matched and passes_size
                truth_matched += int(event_truth)
                any_after_size += int(event_any_size)
                truth_after_size += int(event_truth_size)
            result[threshold] = {
                "signal_any_candidate": any_candidate,
                "signal_truth_matched": truth_matched,
                "signal_any_after_size": any_after_size,
                "signal_truth_after_size": truth_after_size,
            }
    return len(selected), result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--background-volume-glob",
        default="/Users/fabioali/cernbox/CNN/background_scan_volumes/background_scan_*.h5",
    )
    parser.add_argument(
        "--background-prediction-glob",
        default=("runs/cnn21d_signal_t5_center40/background_scan_full_t098/"
                 "background_scan_*_t098.h5"),
    )
    parser.add_argument("--signal-volumes", type=Path, default=Path("/private/tmp/crop_bw_scan_splits/signal_scan_test.h5"))
    parser.add_argument("--signal-predictions", type=Path, default=Path("runs/cnn21d_signal_t5_center40/scan_volumes_test_gt10_t098.h5"))
    parser.add_argument("--thresholds", type=float, nargs="+", default=[0.90, 0.95, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 0.997, 0.999])
    parser.add_argument("--theta-min", type=float, default=10.0)
    parser.add_argument("--min-component-voxels", type=int, default=8)
    parser.add_argument("--truth-propagation-limit", type=int)
    parser.add_argument("--output-dir", type=Path, default=Path("runs/cnn21d_signal_t5_center40/score_size_operating_points"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = sorted(set(args.thresholds))

    background = evaluate_background(
        args.background_volume_glob, args.background_prediction_glob,
        thresholds, args.min_component_voxels,
    )
    signal_total, signal = evaluate_signal(
        args.signal_volumes, args.signal_predictions, thresholds,
        args.theta_min, args.min_component_voxels, args.truth_propagation_limit,
    )
    rows: list[dict[str, object]] = []
    for threshold in thresholds:
        row: dict[str, object] = {"threshold": threshold, "signal_total": signal_total}
        row.update(signal[threshold])
        row.update(background[threshold])
        row["signal_truth_efficiency"] = signal[threshold]["signal_truth_matched"] / signal_total
        row["signal_truth_efficiency_after_size"] = signal[threshold]["signal_truth_after_size"] / signal_total
        row["signal_missed_truth"] = signal_total - signal[threshold]["signal_truth_matched"]
        rows.append(row)

    csv_path = args.output_dir / "score_size_operating_points.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "signal_selection": f"theta_true > {args.theta_min:g} mrad",
        "signal_events": signal_total,
        "size_postfilter": f"largest_component_voxels >= {args.min_component_voxels}",
        "truth_propagation_limit": args.truth_propagation_limit,
        "background_maps": 323,
        "background_windows": 116603,
        "rows": rows,
        "csv": str(csv_path),
    }
    json_path = args.output_dir / "summary.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"csv": str(csv_path.resolve()), "summary": str(json_path.resolve())}, indent=2))


if __name__ == "__main__":
    main()
