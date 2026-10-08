#!/usr/bin/env python3
"""Estimate CNN and Poisson-hot fit angles in original 57x32x32 ROOT crops."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import torch
import uproot

from visualization.make_background_candidate_gifs import estimate_shower_start
from training.train_cnn21d import build_model, load_config


Z_TO_XY_RATIO = 27.0  # 1350 um per XYPseg layer / 50 um per xy bin.


def sliding_positions(length: int, crop_size: int, stride: int) -> list[int]:
    if crop_size <= 0 or stride <= 0 or crop_size > length:
        raise ValueError("invalid crop size or stride")
    positions = list(range(0, length - crop_size + 1, stride))
    if positions[-1] != length - crop_size:
        positions.append(length - crop_size)
    return positions


def theta_mrad(slope_xy: np.ndarray) -> np.ndarray:
    slope = np.asarray(slope_xy, dtype=np.float64)
    return 1000.0 * np.arctan(np.linalg.norm(slope, axis=-1) / Z_TO_XY_RATIO)


def resolve_device(requested: str) -> torch.device:
    if requested == "auto":
        if torch.backends.mps.is_available():
            return torch.device("mps")
        if torch.cuda.is_available():
            return torch.device("cuda")
        return torch.device("cpu")
    return torch.device(requested)


def load_model(config_path: Path, checkpoint_path: Path, device: torch.device) -> torch.nn.Module:
    model = build_model(load_config(config_path)).to(device)
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    return model


def scan_crops(
    raw: np.ndarray,
    background_mu: np.ndarray,
    classifier: torch.nn.Module,
    regression: torch.nn.Module,
    *,
    crop_size: int,
    stride: int,
    batch_size: int,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int]]]:
    """Return score/slope grids flattened over all windows of each 32x32 crop."""
    if raw.ndim != 4 or raw.shape[1:] != (57, 32, 32):
        raise ValueError(f"expected [batch,57,32,32], got {raw.shape}")
    x_positions = sliding_positions(32, crop_size, stride)
    y_positions = sliding_positions(32, crop_size, stride)
    locations = [(x, y) for y in y_positions for x in x_positions]
    crops = np.stack(
        [raw[:, :, y:y + crop_size, x:x + crop_size] for x, y in locations], axis=1
    ).astype(np.float32, copy=False)
    mu = np.asarray(background_mu, dtype=np.float32).reshape(-1, 1, 1, 1, 1)
    crops = (crops - mu) / np.sqrt(mu + np.float32(1e-6))
    flattened = crops.reshape(-1, 57, crop_size, crop_size)
    scores: list[np.ndarray] = []
    slopes: list[np.ndarray] = []
    with torch.no_grad():
        for first in range(0, len(flattened), batch_size):
            batch = torch.from_numpy(flattened[first:first + batch_size, None]).to(device)
            classification = classifier(batch)
            regression_prediction = classification if regression is classifier else regression(batch)
            scores.append(torch.sigmoid(classification["presence_logit"]).cpu().numpy())
            slopes.append(regression_prediction["slope_xy"].cpu().numpy())
    window_count = len(locations)
    return (
        np.concatenate(scores).reshape(len(raw), window_count),
        np.concatenate(slopes).reshape(len(raw), window_count, 2),
        locations,
    )


def fit_record(
    raw: np.ndarray,
    background_mu: float,
    cnn_slope: np.ndarray,
    crop_center_xy: tuple[float, float],
) -> dict[str, object]:
    coordinates = np.arange(raw.shape[-1], dtype=np.float64)
    estimate = estimate_shower_start(
        raw,
        background_mu,
        coordinates,
        coordinates,
        predicted_slope_xy=cnn_slope,
        candidate_center_xy_um=np.asarray(crop_center_xy, dtype=np.float64),
    )
    slope = np.asarray([
        estimate.get("selected_component_slope_x", np.nan),
        estimate.get("selected_component_slope_y", np.nan),
    ], dtype=np.float64)
    has_slope = bool(np.all(np.isfinite(slope)))
    def csv_value(value: object) -> object:
        return "" if value is None else value

    return {
        # A q-hot component is optional.  Keep the row in the output when none
        # is present so that CNN predictions can still be studied on all crops.
        "fit_slope_x_bins_per_layer": float(slope[0]) if has_slope else "",
        "fit_slope_y_bins_per_layer": float(slope[1]) if has_slope else "",
        "fit_theta_mrad": float(theta_mrad(slope)) if has_slope else "",
        "cnn_fit_slope_disagreement_mrad": (
            float(np.linalg.norm(cnn_slope - slope) * 1000.0 / Z_TO_XY_RATIO)
            if has_slope else ""
        ),
        "fit_component_selection": estimate.get("component_selection", "not_found"),
        "fit_start_layer": csv_value(estimate.get("start_layer")),
        "fit_start_x_bin": csv_value(estimate.get("start_x_um")),
        "fit_start_y_bin": csv_value(estimate.get("start_y_um")),
        "fit_quality": estimate.get("start_quality", "not_found"),
        "fit_hot_threshold": estimate.get("poisson_hot_count_threshold", ""),
        "fit_component_voxels": estimate.get("start_component_voxels", 0),
        "fit_component_layers": estimate.get("start_component_layers", 0),
    }


def process_file(
    input_path: Path,
    output_path: Path,
    classifier: torch.nn.Module,
    regression: torch.nn.Module,
    *,
    crop_size: int,
    stride: int,
    batch_size: int,
    chunk_size: int,
    threshold: float,
    max_entries: int | None,
    device: torch.device,
) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    with uproot.open(input_path) as root_file:
        tree = root_file["samples"]
        entries = int(tree.num_entries)
        if max_entries is not None:
            entries = min(entries, max_entries)
        for first in range(0, entries, chunk_size):
            stop = min(entries, first + chunk_size)
            arrays = tree.arrays(
                ["counts", "background_mu", "sample_type", "cell_id", "tag_cell_id"],
                entry_start=first,
                entry_stop=stop,
                library="np",
            )
            raw = np.asarray(arrays["counts"], dtype=np.float32)
            background_mu = np.asarray(arrays["background_mu"], dtype=np.float32)
            scores, slopes, locations = scan_crops(
                raw, background_mu, classifier, regression,
                crop_size=crop_size, stride=stride, batch_size=batch_size, device=device,
            )
            for local_index in range(len(raw)):
                window_index = int(np.argmax(scores[local_index]))
                cnn_slope = np.asarray(slopes[local_index, window_index], dtype=np.float64)
                x_start, y_start = locations[window_index]
                score = float(scores[local_index, window_index])
                row: dict[str, object] = {
                    "root_entry": first + local_index,
                    "sample_type": int(arrays["sample_type"][local_index]),
                    "cell_id": int(arrays["cell_id"][local_index]),
                    "tag_cell_id": int(arrays["tag_cell_id"][local_index]),
                    "background_mu": float(background_mu[local_index]),
                    "cnn_presence_score": score,
                    "cnn_window_x_start": x_start,
                    "cnn_window_y_start": y_start,
                    "cnn_slope_x_bins_per_layer": float(cnn_slope[0]),
                    "cnn_slope_y_bins_per_layer": float(cnn_slope[1]),
                    "cnn_theta_mrad": float(theta_mrad(cnn_slope)),
                    "selected_for_fit": score >= threshold,
                }
                if score >= threshold:
                    row.update(fit_record(
                        raw[local_index], float(background_mu[local_index]), cnn_slope,
                        (x_start + crop_size / 2.0, y_start + crop_size / 2.0),
                    ))
                else:
                    row.update({
                        "fit_slope_x_bins_per_layer": "",
                        "fit_slope_y_bins_per_layer": "",
                        "fit_theta_mrad": "",
                        "cnn_fit_slope_disagreement_mrad": "",
                        "fit_component_selection": "",
                        "fit_start_layer": "",
                        "fit_start_x_bin": "",
                        "fit_start_y_bin": "",
                        "fit_quality": "",
                        "fit_hot_threshold": "",
                        "fit_component_voxels": "",
                        "fit_component_layers": "",
                    })
                rows.append(row)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    selected = [row for row in rows if bool(row["selected_for_fit"])]
    return {
        "input": str(input_path.resolve()),
        "output": str(output_path.resolve()),
        "entries": len(rows),
        "cnn_threshold": threshold,
        "cnn_selected": len(selected),
        "fit_coherent": sum(row.get("fit_quality") == "coherent" for row in selected),
        "cnn_theta_mrad_quantiles_selected": (
            {str(q): float(np.quantile([float(row["cnn_theta_mrad"]) for row in selected], q))
             for q in (0.0, 0.5, 0.9, 1.0)} if selected else {}
        ),
        "fit_theta_mrad_quantiles_coherent": (
            {str(q): float(np.quantile(
                [float(row["fit_theta_mrad"]) for row in selected if row.get("fit_quality") == "coherent"], q
            )) for q in (0.0, 0.5, 0.9, 1.0)}
            if any(row.get("fit_quality") == "coherent" for row in selected) else {}
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--regression-checkpoint", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.90)
    parser.add_argument("--cnn-crop-size", type=int, default=20)
    parser.add_argument("--cnn-stride", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--chunk-size", type=int, default=32)
    parser.add_argument("--max-entries", type=int)
    parser.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="auto")
    args = parser.parse_args()
    if not 0.0 <= args.threshold <= 1.0:
        raise ValueError("threshold must be in [0, 1]")
    if args.max_entries is not None and args.max_entries <= 0:
        raise ValueError("max entries must be positive")

    device = resolve_device(args.device)
    classifier = load_model(args.config, args.checkpoint, device)
    regression = classifier
    if args.regression_checkpoint is not None:
        regression = load_model(args.config, args.regression_checkpoint, device)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for input_path in args.input:
        output_path = args.output_dir / f"{input_path.stem}_angles.csv"
        summaries.append(process_file(
            input_path, output_path, classifier, regression,
            crop_size=args.cnn_crop_size, stride=args.cnn_stride,
            batch_size=args.batch_size, chunk_size=args.chunk_size,
            threshold=args.threshold, max_entries=args.max_entries, device=device,
        ))
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps({
        "config": str(args.config),
        "classification_checkpoint": str(args.checkpoint),
        "regression_checkpoint": str(args.regression_checkpoint or args.checkpoint),
        "cnn_windowing": {
            "source_crop": "57x32x32 ROOT crop",
            "cnn_crop_size": args.cnn_crop_size,
            "cnn_stride": args.cnn_stride,
            "selection": "maximum CNN presence score across windows",
        },
        "fit": "CNN-associated Poisson-hot ridge/component fit in the full 32x32 crop",
        "files": summaries,
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path.resolve()), "files": summaries}, indent=2))


if __name__ == "__main__":
    main()
