#!/usr/bin/env python3
"""Aggregate truth-free background scan shards and export manual-check candidates."""

from __future__ import annotations

import argparse
import csv
import glob
import json
from collections import deque
from pathlib import Path

import h5py
import numpy as np



def connected_components(mask: np.ndarray) -> list[list[tuple[int, int]]]:
    """Return eight-connected components without importing the PyTorch scanner."""
    values = np.asarray(mask, dtype=bool)
    seen = np.zeros_like(values)
    result: list[list[tuple[int, int]]] = []
    for row, col in zip(*np.nonzero(values), strict=True):
        if seen[row, col]:
            continue
        queue = deque([(int(row), int(col))])
        seen[row, col] = True
        component: list[tuple[int, int]] = []
        while queue:
            current_row, current_col = queue.popleft()
            component.append((current_row, current_col))
            for delta_row in (-1, 0, 1):
                for delta_col in (-1, 0, 1):
                    next_row = current_row + delta_row
                    next_col = current_col + delta_col
                    if (
                        0 <= next_row < values.shape[0]
                        and 0 <= next_col < values.shape[1]
                        and values[next_row, next_col]
                        and not seen[next_row, next_col]
                    ):
                        seen[next_row, next_col] = True
                        queue.append((next_row, next_col))
        result.append(component)
    return result


CANDIDATE_FIELDS = (
    "event_id", "run_wall_brick", "brick_name",
    "cell_id", "cell_x", "cell_y", "source_file", "candidate_id",
    "grid_row", "grid_col", "cluster_size", "presence_score",
    "slope_x", "slope_y", "theta_pred_mrad", "phi_pred_rad",
    "x_center_um", "y_center_um",
)


def decode_string(value: object) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-glob", required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument(
        "--thresholds", type=float, nargs="+",
        default=[0.9, 0.95, 0.97, 0.98, 0.985, 0.99, 0.995, 0.997, 0.999],
    )
    parser.add_argument(
        "--exclude-cell-id", type=int, action="append", default=[],
        help="Exclude a cell from all aggregates; repeat for multiple cells",
    )
    args = parser.parse_args()
    excluded_cell_ids = sorted(set(args.exclude_cell_id))

    paths = [Path(path) for path in sorted(glob.glob(args.input_glob))]
    if not paths:
        raise ValueError(f"no prediction files match {args.input_glob!r}")
    operating_points = {
        str(threshold): {"windows": 0, "clusters": 0, "cells": 0}
        for threshold in args.thresholds
    }
    candidates: list[dict[str, object]] = []
    maxima: list[float] = []
    event_ids: list[int] = []
    total_windows = 0
    stored_thresholds: set[float] = set()

    for path in paths:
        with h5py.File(path, "r") as hdf5:
            stored_thresholds.add(float(hdf5.attrs["presence_threshold"]))
            all_scores = np.asarray(hdf5["window_presence_score"])
            shard_cell_ids = np.asarray(hdf5["cell_id"], dtype=np.int64)
            keep = ~np.isin(shard_cell_ids, excluded_cell_ids)
            scores = all_scores[keep]
            maxima.extend(scores.reshape(len(scores), -1).max(axis=1).tolist())
            shard_event_ids = np.asarray(hdf5["event_id"], dtype=np.int64)
            event_ids.extend(shard_event_ids[keep].tolist())
            total_windows += int(scores.size)
            for threshold in args.thresholds:
                cluster_counts = [
                    len(connected_components(score_grid >= threshold))
                    for score_grid in scores
                ]
                point = operating_points[str(threshold)]
                point["windows"] += int((scores >= threshold).sum())
                point["clusters"] += int(sum(cluster_counts))
                point["cells"] += int(sum(count > 0 for count in cluster_counts))

            group = hdf5["candidates"]
            source_files = hdf5["source_file"][:]
            run_wall_bricks = (
                np.asarray(hdf5["run_wall_brick"], dtype=np.int64)
                if "run_wall_brick" in hdf5 else np.full(len(shard_event_ids), -1, dtype=np.int64)
            )
            brick_names = (
                hdf5["brick_name"].asstr()[:]
                if "brick_name" in hdf5 else np.asarray([""] * len(shard_event_ids))
            )
            for index in range(len(group["event_id"])):
                event_index = int(group["event_index"][index])
                if int(shard_cell_ids[event_index]) in excluded_cell_ids:
                    continue
                candidates.append({
                    "event_id": int(group["event_id"][index]),
                    "run_wall_brick": int(run_wall_bricks[event_index]),
                    "brick_name": str(brick_names[event_index]),
                    "cell_id": int(shard_cell_ids[event_index]),
                    "cell_x": int(hdf5["cell_x"][event_index]),
                    "cell_y": int(hdf5["cell_y"][event_index]),
                    "source_file": decode_string(source_files[event_index]),
                    "candidate_id": int(group["candidate_id"][index]),
                    "grid_row": int(group["grid_row"][index]),
                    "grid_col": int(group["grid_col"][index]),
                    "cluster_size": int(group["cluster_size"][index]),
                    "presence_score": float(group["presence_score"][index]),
                    "slope_x": float(group["slope_x"][index]),
                    "slope_y": float(group["slope_y"][index]),
                    "theta_pred_mrad": float(group["theta_pred_mrad"][index]),
                    "phi_pred_rad": float(group["phi_pred_rad"][index]),
                    "x_center_um": float(group["x_center_um"][index]),
                    "y_center_um": float(group["y_center_um"][index]),
                })

    if len(event_ids) != len(set(event_ids)):
        raise ValueError("duplicate event IDs across prediction shards")
    candidates.sort(key=lambda row: float(row["presence_score"]), reverse=True)
    theta = np.asarray([row["theta_pred_mrad"] for row in candidates], dtype=np.float64)
    maxima_array = np.asarray(maxima)
    result = {
        "input_glob": args.input_glob,
        "excluded_cell_ids": excluded_cell_ids,
        "prediction_files": len(paths),
        "cells": len(event_ids),
        "windows": total_windows,
        "stored_candidate_thresholds": sorted(stored_thresholds),
        "stored_candidates": len(candidates),
        "cells_with_stored_candidate": len({int(row["event_id"]) for row in candidates}),
        "operating_points": operating_points,
        "cell_maximum_score_quantiles": {
            str(q): float(np.quantile(maxima_array, q))
            for q in (0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
        },
        "stored_candidate_theta_bins": {
            "lt_10": int((theta < 10).sum()),
            "10_to_20": int(((theta >= 10) & (theta < 20)).sum()),
            "20_to_50": int(((theta >= 20) & (theta < 50)).sum()),
            "ge_50": int((theta >= 50).sum()),
        },
    }

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_FIELDS)
        writer.writeheader()
        writer.writerows(candidates)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
