#!/usr/bin/env python3
"""Compare manually identified data showers with every CNN scan window."""
from __future__ import annotations

import argparse
import csv
import glob
import re
from pathlib import Path

import h5py
import numpy as np


def decode_direction(code: int) -> tuple[float, float]:
    """Decode ((-ty+50)/2)*51 + ((-tx+50)/2)."""
    row, col = divmod(code, 51)
    return 50.0 - 2.0 * col, 50.0 - 2.0 * row


def parse_rows(path: Path) -> list[dict[str, float | int]]:
    rows = []
    for line_no, line in enumerate(path.read_text().splitlines(), start=1):
        fields = [item.strip() for item in line.split("*") if item.strip()]
        if not fields or line.lstrip().startswith("#"):
            continue
        if len(fields) != 4:
            raise ValueError(f"{path}:{line_no}: expected code * cell * x_um * y_um")
        code, cell, x, y = int(fields[0]), int(fields[1]), float(fields[2]), float(fields[3])
        tx, ty = decode_direction(code)
        rows.append({"manual_code": code, "cell_id": cell, "x_um": x, "y_um": y,
                     "manual_tx_mrad": tx, "manual_ty_mrad": ty,
                     "manual_theta_mrad": float(np.hypot(tx, ty))})
    return rows


def source_for_cell(volume_glob: str, cell_id: int) -> Path:
    for path in map(Path, glob.glob(volume_glob)):
        if path.suffix != ".h5":
            continue
        with h5py.File(path, "r") as f:
            if cell_id in set(np.asarray(f["cell_id"], dtype=int)):
                return path
    raise ValueError(f"no input volume contains cell {cell_id}")


def prediction_for_volume(prediction_glob: str, volume: Path) -> Path:
    prefix = volume.stem.split("_cell", 1)[0]
    candidates = [Path(p) for p in glob.glob(prediction_glob)
                  if Path(p).name.startswith(prefix)]
    if len(candidates) != 1:
        raise ValueError(f"cannot uniquely match prediction to {volume.name}: {candidates}")
    return candidates[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manual-list", type=Path, required=True)
    parser.add_argument("--volume-glob", required=True)
    parser.add_argument("--prediction-glob", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--local-radius-um", type=float, default=1500.0)
    args = parser.parse_args()

    output = []
    for item in parse_rows(args.manual_list):
        volume = source_for_cell(args.volume_glob, int(item["cell_id"]))
        pred_path = prediction_for_volume(args.prediction_glob, volume)
        with h5py.File(volume, "r") as source, h5py.File(pred_path, "r") as pred:
            indices = np.flatnonzero(np.asarray(source["cell_id"], dtype=int) == item["cell_id"])
            if len(indices) != 1:
                raise ValueError(f"expected one event for cell {item['cell_id']}")
            event_index = int(indices[0])
            x_edges = np.asarray(source["x_edges_um"][event_index])
            y_edges = np.asarray(source["y_edges_um"][event_index])
            x_grid = np.asarray(pred["grid_x_start_bin"], dtype=int)
            y_grid = np.asarray(pred["grid_y_start_bin"], dtype=int)
            crop_size = int(pred.attrs["crop_size"])
            x_center = 0.5 * (x_edges[x_grid] + x_edges[x_grid + crop_size])
            y_center = 0.5 * (y_edges[y_grid] + y_edges[y_grid + crop_size])
            xx, yy = np.meshgrid(x_center, y_center)
            distance = np.hypot(xx - item["x_um"], yy - item["y_um"])
            score = np.asarray(pred["window_presence_score"][event_index])
            slope = np.asarray(pred["window_slope_xy"][event_index])
            nearest = np.unravel_index(np.argmin(distance), distance.shape)
            local = distance <= args.local_radius_um
            best = np.unravel_index(np.argmax(np.where(local, score, -np.inf)), score.shape)
            global_best = np.unravel_index(np.argmax(score), score.shape)
            for key, index in (("nearest", nearest), ("local_best", best), ("global_best", global_best)):
                row, col = index
                item[f"{key}_score"] = float(score[row, col])
                item[f"{key}_distance_um"] = float(distance[row, col])
                item[f"{key}_grid_row"] = int(row)
                item[f"{key}_grid_col"] = int(col)
                item[f"{key}_tx_mrad"] = float(1000.0 * slope[row, col, 0] / 27.0)
                item[f"{key}_ty_mrad"] = float(1000.0 * slope[row, col, 1] / 27.0)
            item["volume"] = volume.name
            item["prediction"] = pred_path.name
            item["distance_to_x_edge_um"] = float(min(item["x_um"] - x_edges[0], x_edges[-1] - item["x_um"]))
            item["distance_to_y_edge_um"] = float(min(item["y_um"] - y_edges[0], y_edges[-1] - item["y_um"]))
            output.append(item)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(output[0]))
        writer.writeheader(); writer.writerows(output)
    found = sum(float(row["local_best_score"]) >= 0.9 for row in output)
    print(f"wrote {args.output}: {found}/{len(output)} local matches have score >= 0.9")


if __name__ == "__main__":
    main()
