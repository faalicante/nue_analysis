#!/usr/bin/env python3
"""Export target scan candidates that do not overlap a reference model."""

from __future__ import annotations

import argparse
import csv
import glob
import json
from pathlib import Path

import h5py
import numpy as np

from candidates.export_background_scan_candidates import CANDIDATE_FIELDS, connected_components, decode_string


def components_by_event(paths: list[Path], threshold: float) -> dict[int, list[set[tuple[int, int]]]]:
    """Map every event to its eight-connected score components."""
    result: dict[int, list[set[tuple[int, int]]]] = {}
    for path in paths:
        with h5py.File(path, "r") as hdf5:
            for event_id, scores in zip(
                hdf5["event_id"][:], hdf5["window_presence_score"][:], strict=True,
            ):
                key = int(event_id)
                if key in result:
                    raise ValueError(f"duplicate event ID {key} in {path}")
                result[key] = [set(component) for component in connected_components(scores >= threshold)]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-prediction-glob", required=True)
    parser.add_argument("--reference-prediction-glob", required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.90)
    args = parser.parse_args()
    if not 0.0 <= args.threshold <= 1.0:
        parser.error("--threshold must be in [0, 1]")

    target_paths = sorted(map(Path, glob.glob(args.target_prediction_glob)))
    reference_paths = sorted(map(Path, glob.glob(args.reference_prediction_glob)))
    if not target_paths or not reference_paths:
        parser.error("target or reference prediction glob is empty")

    reference = components_by_event(reference_paths, args.threshold)
    candidates: list[dict[str, object]] = []
    total_target_candidates = 0
    overlap_excluded = 0

    for path in target_paths:
        with h5py.File(path, "r") as hdf5:
            group = hdf5["candidates"]
            scores = hdf5["window_presence_score"]
            cell_ids = np.asarray(hdf5["cell_id"], dtype=np.int64)
            cell_x = np.asarray(hdf5["cell_x"], dtype=np.int64)
            cell_y = np.asarray(hdf5["cell_y"], dtype=np.int64)
            source_files = hdf5["source_file"][:]
            event_ids = np.asarray(hdf5["event_id"], dtype=np.int64)
            brick_names = (
                hdf5["brick_name"].asstr()[:]
                if "brick_name" in hdf5 else np.asarray([""] * len(event_ids))
            )
            run_wall_bricks = (
                np.asarray(hdf5["run_wall_brick"], dtype=np.int64)
                if "run_wall_brick" in hdf5 else np.full(len(event_ids), -1, dtype=np.int64)
            )
            for index in range(len(group["event_id"])):
                total_target_candidates += 1
                event_index = int(group["event_index"][index])
                row = int(group["grid_row"][index])
                col = int(group["grid_col"][index])
                component = next(
                    (
                        set(points)
                        for points in connected_components(scores[event_index] >= args.threshold)
                        if (row, col) in points
                    ),
                    set(),
                )
                event_id = int(event_ids[event_index])
                if any(component & other for other in reference.get(event_id, [])):
                    overlap_excluded += 1
                    continue
                candidates.append({
                    "event_id": event_id,
                    "run_wall_brick": int(run_wall_bricks[event_index]),
                    "brick_name": str(brick_names[event_index]),
                    "cell_id": int(cell_ids[event_index]),
                    "cell_x": int(cell_x[event_index]),
                    "cell_y": int(cell_y[event_index]),
                    "source_file": decode_string(source_files[event_index]),
                    "candidate_id": int(group["candidate_id"][index]),
                    "grid_row": row,
                    "grid_col": col,
                    "cluster_size": int(group["cluster_size"][index]),
                    "presence_score": float(group["presence_score"][index]),
                    "slope_x": float(group["slope_x"][index]),
                    "slope_y": float(group["slope_y"][index]),
                    "theta_pred_mrad": float(group["theta_pred_mrad"][index]),
                    "phi_pred_rad": float(group["phi_pred_rad"][index]),
                    "x_center_um": float(group["x_center_um"][index]),
                    "y_center_um": float(group["y_center_um"][index]),
                })

    candidates.sort(key=lambda row: float(row["presence_score"]), reverse=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_FIELDS)
        writer.writeheader()
        writer.writerows(candidates)
    payload = {
        "threshold": args.threshold,
        "target_prediction_glob": args.target_prediction_glob,
        "reference_prediction_glob": args.reference_prediction_glob,
        "target_candidate_count": total_target_candidates,
        "overlap_excluded_candidate_count": overlap_excluded,
        "exclusive_candidate_count": len(candidates),
        "cells_with_exclusive_candidate": len({int(row["event_id"]) for row in candidates}),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
