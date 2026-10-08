#!/usr/bin/env python3
"""Fit straight shower axes directly to significant voxels in the raw HDF5 volumes."""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/cnn21d-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
import numpy as np

from training.cnn_dataset import CNNDataConfig, HDF5CNN21DDataset
from training.train_cnn21d import _stratified_local_indices, load_config


@dataclass
class LineFit:
    slope_x: float
    slope_y: float
    intercept_x: float
    intercept_y: float
    theta_mrad: float
    points: int
    slices: int
    z_span: int
    radial_rms: float
    score: float
    threshold_points: int
    reliable: bool


def fit_axis_ransac(
    raw: np.ndarray,
    background_mu: np.ndarray,
    *,
    sigma_threshold: float = 5.0,
    radius: float = 1.6,
    iterations: int = 10_000,
    seed: int = 1,
) -> tuple[LineFit | None, np.ndarray]:
    mu = np.asarray(background_mu, dtype=np.float64)
    significance = (raw.astype(np.float64) - mu[:, None, None]) / np.sqrt(mu[:, None, None] + 1e-6)
    z, y, x = np.where(significance > sigma_threshold)
    if len(z) < 2:
        return None, significance
    weights = np.clip(significance[z, y, x] - sigma_threshold + 0.25, 0.25, 5.0)
    rng = np.random.default_rng(seed)
    best: tuple[float, np.ndarray] | None = None
    for _ in range(iterations):
        first, second = rng.integers(0, len(z), size=2)
        delta_z = int(z[second]) - int(z[first])
        if abs(delta_z) < 5:
            continue
        slope_x = (x[second] - x[first]) / delta_z
        slope_y = (y[second] - y[first]) / delta_z
        magnitude = float(np.hypot(slope_x, slope_y))
        if not 0.05 < magnitude < 2.0:
            continue
        predicted_x = x[first] + slope_x * (z - z[first])
        predicted_y = y[first] + slope_y * (z - z[first])
        mask = np.hypot(x - predicted_x, y - predicted_y) < radius
        unique_slices = np.unique(z[mask]).size
        span = int(np.ptp(z[mask])) if mask.any() else 0
        score = float(weights[mask].sum() + 1.5 * unique_slices + 0.1 * span)
        if best is None or score > best[0]:
            best = score, mask
    if best is None:
        return None, significance
    score, mask = best
    selected_z = z[mask].astype(np.float64)
    design = np.column_stack((np.ones(mask.sum()), selected_z))
    selected_weights = weights[mask]
    weighted_design = design * np.sqrt(selected_weights)[:, None]
    fit_x = np.linalg.lstsq(weighted_design, x[mask] * np.sqrt(selected_weights), rcond=None)[0]
    fit_y = np.linalg.lstsq(weighted_design, y[mask] * np.sqrt(selected_weights), rcond=None)[0]
    predicted_x = design @ fit_x; predicted_y = design @ fit_y
    rms = float(np.sqrt(np.average(
        (x[mask] - predicted_x) ** 2 + (y[mask] - predicted_y) ** 2,
        weights=selected_weights,
    )))
    slices = int(np.unique(z[mask]).size); z_span = int(np.ptp(z[mask]))
    slope_x, slope_y = float(fit_x[1]), float(fit_y[1])
    theta = float(1000.0 * np.arctan(np.hypot(slope_x, slope_y) / 27.0))
    result = LineFit(
        slope_x=slope_x, slope_y=slope_y,
        intercept_x=float(fit_x[0]), intercept_y=float(fit_y[0]), theta_mrad=theta,
        points=int(mask.sum()), slices=slices, z_span=z_span, radial_rms=rms,
        score=float(score), threshold_points=len(z),
        reliable=slices >= 8 and z_span >= 8 and rms <= 1.25,
    )
    return result, significance


def selected_pilot_indices(config: dict[str, Any], split: str) -> np.ndarray:
    data_config = CNNDataConfig(**config["data"]["dataset"], seed=int(config["seed"]))
    dataset = HDF5CNN21DDataset(config["data"]["hdf5_path"], split=split, training=False, config=data_config)  # type: ignore[arg-type]
    task = config["task"]
    theta_range = tuple(float(x) for x in task["signal_theta_range_mrad"])
    local = _stratified_local_indices(
        dataset, config["subsets"]["pilot"][split],
        int(config["seed"]) + {"train": 0, "validation": 1, "test": 2}[split],
        tuple(int(x) for x in task["include_sample_types"]),
        tuple(int(x) for x in task["positive_sample_types"]),
        theta_range,
    )
    return dataset.indices[np.asarray(local, dtype=np.int64)]


def save_fit_plot(
    path: Path,
    hdf5_index: int,
    signal_id: int,
    label: np.ndarray,
    fit: LineFit | None,
    significance: np.ndarray,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    z = np.arange(significance.shape[0])
    projections = (significance.max(axis=1).T, significance.max(axis=2).T)
    names = ("x", "y")
    coordinate_sizes = (significance.shape[2], significance.shape[1])
    label_slopes = (float(label[0]), float(label[1]))
    if fit is not None:
        fit_slopes = (fit.slope_x, fit.slope_y)
        fit_intercepts = (fit.intercept_x, fit.intercept_y)
        # Place the target line through the weighted center of the fitted line;
        # only its slope, not this nuisance intercept, is compared.
        center_z = 0.5 * (significance.shape[0] - 1)
        target_intercepts = tuple(
            fit_intercepts[index] + (fit_slopes[index] - label_slopes[index]) * center_z
            for index in range(2)
        )
    for index, (axis, projection, name) in enumerate(zip(axes, projections, names, strict=True)):
        axis.imshow(projection, origin="lower", aspect="auto", cmap="magma", extent=(0, 56, 0, coordinate_sizes[index]))
        if fit is not None:
            axis.plot(z, fit_intercepts[index] + fit_slopes[index] * z, color="lime", linewidth=2, label="fit libero")
            axis.plot(z, target_intercepts[index] + label_slopes[index] * z, color="cyan", linewidth=2, linestyle="--", label="slope label")
        axis.set(xlabel="slice z", ylabel=f"coordinata {name} [bin]", title=f"proiezione z–{name}")
        axis.set_xlim(0, 56); axis.legend(fontsize=8)
    fit_text = "fit non disponibile" if fit is None else (
        f"fit=({fit.slope_x:.3f},{fit.slope_y:.3f}), theta={fit.theta_mrad:.1f} mrad, "
        f"slice={fit.slices}, RMS={fit.radial_rms:.2f}"
    )
    theta_label = 1000.0 * np.arctan(np.linalg.norm(label) / 27.0)
    figure.suptitle(
        f"HDF5 {hdf5_index} · signal {signal_id} · label=({label[0]:.3f},{label[1]:.3f}), theta={theta_label:.1f} mrad\n{fit_text}"
    )
    figure.savefig(path, dpi=170); plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("training/configs/cnn21d.yaml"))
    parser.add_argument("--output-dir", type=Path, default=Path("runs/cnn21d_sample_type_no_poisson_theta10_50/pilot/centroid_fit"))
    parser.add_argument("--iterations", type=int, default=10_000)
    args = parser.parse_args(); config = load_config(args.config)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    indices = np.concatenate((selected_pilot_indices(config, "validation"), selected_pilot_indices(config, "test")))
    selected_examples = {12130, 8335, 9506, 9201, 9876, 12042}
    results: list[dict[str, Any]] = []
    with h5py.File(config["data"]["hdf5_path"], "r") as hdf5:
        for raw_index in indices:
            if int(hdf5["sample_type"][raw_index]) != 2:
                continue
            raw = np.asarray(hdf5["regions_raw"][raw_index])
            background = np.asarray(hdf5["background_mu"][raw_index])
            label = np.asarray(hdf5["slope_xy"][raw_index], dtype=np.float64)
            fit, significance = fit_axis_ransac(
                raw, background, iterations=args.iterations, seed=int(raw_index)
            )
            theta_label = float(1000.0 * np.arctan(np.linalg.norm(label) / 27.0))
            record: dict[str, Any] = {
                "hdf5_index": int(raw_index),
                "signal_event_id": int(hdf5["signal_event_id"][raw_index]),
                "label_slope_xy": label.tolist(),
                "label_theta_mrad": theta_label,
                "raw32_fit": None if fit is None else asdict(fit),
            }
            if fit is not None:
                record["slope_l2_error"] = float(np.linalg.norm(np.asarray([fit.slope_x, fit.slope_y]) - label))
                record["theta_abs_error_mrad"] = abs(fit.theta_mrad - theta_label)
            if int(raw_index) in selected_examples:
                crop_fit, _ = fit_axis_ransac(
                    raw[:, 6:26, 6:26], background,
                    iterations=args.iterations, seed=int(raw_index),
                )
                record["crop20_fit"] = None if crop_fit is None else asdict(crop_fit)
                save_fit_plot(
                    args.output_dir / f"hdf5_{int(raw_index)}_signal_{record['signal_event_id']}.png",
                    int(raw_index), record["signal_event_id"], label, fit, significance,
                )
            results.append(record)
    reliable = [record for record in results if record["raw32_fit"] and record["raw32_fit"]["reliable"]]
    summary = {
        "method": {
            "volume": "raw 32x32",
            "voxel_threshold_sigma": 5.0,
            "ransac_corridor_radius_bins": 1.6,
            "reliable_definition": ">=8 slices, z span >=8, radial RMS <=1.25 bins",
        },
        "signal_count": len(results),
        "reliable_fit_count": len(reliable),
        "reliable_fit_fraction": len(reliable) / len(results) if results else 0.0,
        "reliable_median_slope_l2_error": float(np.median([record["slope_l2_error"] for record in reliable])) if reliable else None,
        "reliable_median_theta_abs_error_mrad": float(np.median([record["theta_abs_error_mrad"] for record in reliable])) if reliable else None,
        "records": results,
    }
    with (args.output_dir / "fit_results.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2)
    print(json.dumps({key: value for key, value in summary.items() if key != "records"}, indent=2))


if __name__ == "__main__":
    main()
