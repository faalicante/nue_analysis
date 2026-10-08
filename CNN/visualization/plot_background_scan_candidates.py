#!/usr/bin/env python3
"""Render 57-layer contact sheets for truth-free background scan candidates."""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from training.cnn_dataset import normalize_per_slice


def candidate_filename(rank: int, row: dict[str, object]) -> str:
    score = f"{float(row['presence_score']):.6f}".replace(".", "p")
    return (
        f"rank_{rank:03d}_cell{int(row['cell_id']):03d}_"
        f"x{int(row['cell_x']):02d}_y{int(row['cell_y']):02d}_"
        f"cand{int(row['candidate_id']):02d}_score{score}.png"
    )


def add_slope_arrow(axis: plt.Axes, slope_xy: np.ndarray) -> float:
    norm = float(np.linalg.norm(slope_xy))
    scale = min(10.0, 6.5 / norm) if norm > 0 else 10.0
    delta = slope_xy * scale
    axis.arrow(
        9.5, 9.5, float(delta[0]), float(delta[1]),
        color="lime", width=0.10, head_width=0.70,
        length_includes_head=True, zorder=3,
    )
    return scale


def save_sheet(record: dict[str, object], output: Path) -> None:
    volume = np.asarray(record["volume_normalized"])
    positive = np.maximum(volume, 0.0)
    vmax = max(float(np.quantile(positive, 0.995)), 1.0)
    vmin = max(float(np.quantile(volume, 0.02)), -3.0)
    slope = np.asarray(record["slope_xy"], dtype=np.float64)
    figure, axes = plt.subplots(8, 8, figsize=(12, 12), constrained_layout=True)
    arrow_scale = 10.0
    for index, axis in enumerate(axes.flat):
        axis.set_xticks([])
        axis.set_yticks([])
        if index < 57:
            axis.imshow(
                volume[index], origin="lower", cmap="magma",
                vmin=vmin, vmax=vmax, interpolation="nearest",
            )
            arrow_scale = add_slope_arrow(axis, slope)
            axis.set_title(f"z={index + 1:02d}", fontsize=7)
        elif index == 57:
            projection = positive.sum(axis=0)
            projection_max = max(float(np.quantile(projection, 0.995)), 1.0)
            axis.imshow(
                projection, origin="lower", cmap="magma",
                vmin=0.0, vmax=projection_max, interpolation="nearest",
            )
            add_slope_arrow(axis, slope)
            axis.set_title("somma z", fontsize=7)
        else:
            axis.axis("off")

    figure.suptitle(
        f"Rank {int(record['rank']):02d} | cell_id={int(record['cell_id'])} "
        f"(x={int(record['cell_x'])}, y={int(record['cell_y'])}) | "
        f"candidate={int(record['candidate_id'])} | "
        f"P(signal)={float(record['presence_score']):.6f} | "
        f"theta_pred={float(record['theta_pred_mrad']):.2f} mrad | "
        f"phi_pred={float(record['phi_pred_rad']):.2f} rad\n"
        f"crop CNN normalizzato 20x20 | freccia verde: slope predetta x{arrow_scale:.1f} layer",
        fontsize=11,
    )
    figure.savefig(output, dpi=150)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--volume-glob", required=True)
    parser.add_argument("--prediction-glob", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    volume_paths = {path.name: path for path in map(Path, glob.glob(args.volume_glob))}
    prediction_paths = sorted(map(Path, glob.glob(args.prediction_glob)))
    if not volume_paths or not prediction_paths:
        raise ValueError("volume or prediction glob is empty")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []

    for prediction_path in prediction_paths:
        volume_name = prediction_path.name.removesuffix("_t098.h5") + ".h5"
        if volume_name not in volume_paths:
            raise ValueError(f"cannot match {prediction_path} to a source volume shard")
        volume_path = volume_paths[volume_name]
        with h5py.File(volume_path, "r") as source, h5py.File(prediction_path, "r") as prediction:
            group = prediction["candidates"]
            x_positions = np.asarray(prediction["grid_x_start_bin"])
            y_positions = np.asarray(prediction["grid_y_start_bin"])
            crop_size = int(prediction.attrs["crop_size"])
            for candidate_index in range(len(group["event_id"])):
                event_index = int(group["event_index"][candidate_index])
                grid_row = int(group["grid_row"][candidate_index])
                grid_col = int(group["grid_col"][candidate_index])
                x_start = int(x_positions[grid_col])
                y_start = int(y_positions[grid_row])
                raw = np.asarray(source["volumes_raw"][event_index])[
                    :, y_start:y_start + crop_size, x_start:x_start + crop_size
                ]
                normalized = normalize_per_slice(raw, float(source["background_mu"][event_index]))
                records.append({
                    "cell_id": int(group["event_id"][candidate_index]),
                    "cell_x": int(source["cell_x"][event_index]),
                    "cell_y": int(source["cell_y"][event_index]),
                    "candidate_id": int(group["candidate_id"][candidate_index]),
                    "grid_row": grid_row,
                    "grid_col": grid_col,
                    "presence_score": float(group["presence_score"][candidate_index]),
                    "slope_xy": np.asarray(group["slope_x"][candidate_index:candidate_index + 1].tolist() + group["slope_y"][candidate_index:candidate_index + 1].tolist()),
                    "theta_pred_mrad": float(group["theta_pred_mrad"][candidate_index]),
                    "phi_pred_rad": float(group["phi_pred_rad"][candidate_index]),
                    "x_center_um": float(group["x_center_um"][candidate_index]),
                    "y_center_um": float(group["y_center_um"][candidate_index]),
                    "background_mu": float(source["background_mu"][event_index]),
                    "source_hdf5": str(volume_path),
                    "source_root": (
                        source["source_file"][event_index].decode("utf-8")
                        if isinstance(source["source_file"][event_index], bytes)
                        else str(source["source_file"][event_index])
                    ),
                    "volume_normalized": normalized,
                })

    records.sort(key=lambda row: float(row["presence_score"]), reverse=True)
    manifest_records: list[dict[str, object]] = []
    for rank, record in enumerate(records, start=1):
        record["rank"] = rank
        path = args.output_dir / candidate_filename(rank, record)
        save_sheet(record, path)
        manifest_records.append({
            key: (value.tolist() if isinstance(value, np.ndarray) else value)
            for key, value in record.items() if key != "volume_normalized"
        } | {"image": str(path)})

    manifest = {
        "candidate_count": len(records),
        "display": {
            "crop": "20x20 representative CNN window",
            "normalization": "(raw - background_mu) / sqrt(background_mu)",
            "layer_numbering": "1..57",
            "arrow": "predicted slope; scale capped to fit the panel",
        },
        "candidates": manifest_records,
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "candidate_count": len(records),
        "output_dir": str(args.output_dir.resolve()),
        "manifest": str((args.output_dir / "manifest.json").resolve()),
    }, indent=2))


if __name__ == "__main__":
    main()

