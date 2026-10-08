#!/usr/bin/env python3
"""Render truth-centred CNN windows for signal events missed by the stride scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from training.cnn_dataset import normalize_per_slice
from training.cnn21d.metrics import slopes_to_angles
from scanning.scan_cnn21d_volumes import cluster_contains_truth, connected_components


def arrow_scale(true_slope: np.ndarray, predicted_slope: np.ndarray) -> float:
    maximum = max(float(np.linalg.norm(true_slope)), float(np.linalg.norm(predicted_slope)), 1e-9)
    return min(10.0, 6.5 / maximum)


def add_arrows(
    axis: plt.Axes,
    true_slope: np.ndarray,
    predicted_slope: np.ndarray,
    scale: float,
) -> None:
    center = (9.5, 9.5)
    true_delta = true_slope * scale
    predicted_delta = predicted_slope * scale
    axis.arrow(
        *center, float(true_delta[0]), float(true_delta[1]),
        color="cyan", width=0.09, head_width=0.65,
        length_includes_head=True, zorder=4,
    )
    axis.arrow(
        *center, float(predicted_delta[0]), float(predicted_delta[1]),
        color="lime", width=0.09, head_width=0.65,
        length_includes_head=True, zorder=5,
    )


def save_sheet(record: dict[str, object], output: Path) -> None:
    volume = np.asarray(record["volume_normalized"])
    positive = np.maximum(volume, 0.0)
    vmin = max(float(np.quantile(volume, 0.02)), -3.0)
    vmax = max(float(np.quantile(positive, 0.995)), 1.0)
    true_slope = np.asarray(record["true_slope"])
    predicted_slope = np.asarray(record["predicted_slope"])
    scale = arrow_scale(true_slope, predicted_slope)
    truth_local = np.asarray(record["truth_local_xy"])
    figure, axes = plt.subplots(8, 8, figsize=(12, 12), constrained_layout=True)
    for index, axis in enumerate(axes.flat):
        axis.set_xticks([])
        axis.set_yticks([])
        if index < 57:
            axis.imshow(
                volume[index], origin="lower", cmap="magma",
                vmin=vmin, vmax=vmax, interpolation="nearest",
            )
            add_arrows(axis, true_slope, predicted_slope, scale)
            axis.add_patch(Circle(
                (float(truth_local[0]), float(truth_local[1])),
                radius=0.7, fill=False, edgecolor="cyan", linewidth=0.8,
            ))
            axis.set_title(f"z={index + 1:02d}", fontsize=7)
        elif index == 57:
            projection = positive.sum(axis=0)
            projection_max = max(float(np.quantile(projection, 0.995)), 1.0)
            axis.imshow(
                projection, origin="lower", cmap="magma",
                vmin=0.0, vmax=projection_max, interpolation="nearest",
            )
            add_arrows(axis, true_slope, predicted_slope, scale)
            axis.add_patch(Circle(
                (float(truth_local[0]), float(truth_local[1])),
                radius=0.8, fill=False, edgecolor="cyan", linewidth=1.0,
            ))
            axis.set_title("somma z", fontsize=7)
        else:
            axis.axis("off")

    off_truth = " | candidato alto fuori truth" if bool(record["has_off_truth_candidate"]) else ""
    figure.suptitle(
        f"Signal non rilevato | event_id={int(record['event_id'])} | "
        f"theta_true={float(record['theta_true_mrad']):.2f} mrad | "
        f"best score su truth={float(record['best_truth_score']):.6f} | "
        f"max score globale={float(record['maximum_score']):.6f}{off_truth}\n"
        f"theta_pred nel best crop={float(record['theta_pred_mrad']):.2f} mrad | "
        f"cyan=vera, verde=predetta (x{scale:.1f} layer), cerchio=centro MC propagato",
        fontsize=10.5,
    )
    figure.savefig(output, dpi=150)
    plt.close(figure)


def save_overview(records: list[dict[str, object]], output: Path) -> None:
    figure, axes = plt.subplots(4, 4, figsize=(13, 13), constrained_layout=True)
    for axis, record in zip(axes.flat, records):
        volume = np.maximum(np.asarray(record["volume_normalized"]), 0.0)
        projection = volume.sum(axis=0)
        vmax = max(float(np.quantile(projection, 0.995)), 1.0)
        true_slope = np.asarray(record["true_slope"])
        predicted_slope = np.asarray(record["predicted_slope"])
        scale = arrow_scale(true_slope, predicted_slope)
        truth_local = np.asarray(record["truth_local_xy"])
        axis.imshow(projection, origin="lower", cmap="magma", vmin=0.0, vmax=vmax)
        add_arrows(axis, true_slope, predicted_slope, scale)
        axis.add_patch(Circle(
            (float(truth_local[0]), float(truth_local[1])),
            radius=0.8, fill=False, edgecolor="cyan", linewidth=1.0,
        ))
        axis.set_title(
            f"evt {int(record['event_id'])} | theta={float(record['theta_true_mrad']):.1f}\n"
            f"score truth={float(record['best_truth_score']):.3f} | max={float(record['maximum_score']):.3f}",
            fontsize=9,
        )
        axis.set_xticks([])
        axis.set_yticks([])
    figure.suptitle(
        "16 signal theta>10 mrad non rilevati — proiezione sui 57 layer\n"
        "cyan=slope vera, verde=predetta, cerchio=centro MC propagato",
        fontsize=12,
    )
    figure.savefig(output, dpi=160)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--volumes", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.98)
    parser.add_argument("--theta-min", type=float, default=10.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, object]] = []
    with h5py.File(args.volumes, "r") as source, h5py.File(args.predictions, "r") as prediction:
        scores = prediction["window_presence_score"][:]
        slopes = prediction["window_slope_xy"][:]
        x_positions = prediction["grid_x_start_bin"][:].astype(int).tolist()
        y_positions = prediction["grid_y_start_bin"][:].astype(int).tolist()
        crop_size = int(prediction.attrs["crop_size"])
        truth_x = source["truth/propagated_center_x_bin"][:]
        truth_y = source["truth/propagated_center_y_bin"][:]
        true_slopes = source["truth/slope_xy_bins_per_z"][:]
        theta_true = slopes_to_angles(true_slopes)[0]

        for event_index in range(len(scores)):
            if theta_true[event_index] <= args.theta_min:
                continue
            components = connected_components(scores[event_index] >= args.threshold)
            detected = any(
                cluster_contains_truth(
                    component,
                    float(truth_x[event_index]), float(truth_y[event_index]),
                    x_positions, y_positions, crop_size,
                )
                for component in components
            )
            if detected:
                continue
            containing = [
                (row, col)
                for row, y0 in enumerate(y_positions)
                for col, x0 in enumerate(x_positions)
                if x0 <= truth_x[event_index] < x0 + crop_size
                and y0 <= truth_y[event_index] < y0 + crop_size
            ]
            if not containing:
                raise ValueError(f"truth for event index {event_index} is outside every scan window")
            row, col = max(containing, key=lambda rc: float(scores[event_index][rc]))
            x0 = x_positions[col]
            y0 = y_positions[row]
            raw = source["volumes_raw"][event_index, :, y0:y0 + crop_size, x0:x0 + crop_size]
            predicted_slope = np.asarray(slopes[event_index, row, col])
            theta_pred = slopes_to_angles(predicted_slope[None])[0][0]
            records.append({
                "event_index": event_index,
                "event_id": int(source["event_id"][event_index]),
                "theta_true_mrad": float(theta_true[event_index]),
                "theta_pred_mrad": float(theta_pred),
                "true_slope": np.asarray(true_slopes[event_index]),
                "predicted_slope": predicted_slope,
                "best_truth_score": float(scores[event_index, row, col]),
                "maximum_score": float(scores[event_index].max()),
                "has_off_truth_candidate": bool(components),
                "grid_row": row,
                "grid_col": col,
                "truth_local_xy": np.asarray([
                    truth_x[event_index] - x0,
                    truth_y[event_index] - y0,
                ]),
                "volume_normalized": normalize_per_slice(
                    raw, float(source["background_mu"][event_index])
                ),
            })

    records.sort(key=lambda record: float(record["theta_true_mrad"]), reverse=True)
    manifest_rows: list[dict[str, object]] = []
    for rank, record in enumerate(records, start=1):
        filename = (
            f"rank_{rank:02d}_event_{int(record['event_id']):05d}_"
            f"theta{float(record['theta_true_mrad']):06.2f}_"
            f"score{float(record['best_truth_score']):.6f}".replace(".", "p")
            + ".png"
        )
        output = args.output_dir / filename
        save_sheet(record, output)
        manifest_rows.append({
            key: value.tolist() if isinstance(value, np.ndarray) else value
            for key, value in record.items() if key != "volume_normalized"
        } | {"image": str(output)})

    overview_path = args.output_dir / "overview_16_missed_signal.png"
    save_overview(records, overview_path)

    manifest = {
        "threshold": args.threshold,
        "theta_min_mrad_strict": args.theta_min,
        "missed_events": len(records),
        "overview": str(overview_path),
        "events": manifest_rows,
    }
    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "missed_events": len(records),
        "output_dir": str(args.output_dir.resolve()),
        "manifest": str(manifest_path.resolve()),
    }, indent=2))


if __name__ == "__main__":
    main()
