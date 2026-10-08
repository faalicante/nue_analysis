#!/usr/bin/env python3
"""Compare angle-threshold pilot operating points and regression by true theta."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from training.cnn_dataset import CNNDataConfig
from training.cnn21d.losses import MultiTaskLoss
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, load_config, make_loader, predict


PILOTS = {
    "slope_weight_0p3": (
        Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20.yaml"),
        Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20/pilot/best_model.pt"),
    ),
    "slope_weight_1p0": (
        Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20_slope1.yaml"),
        Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20_slope1/pilot/best_model.pt"),
    ),
    "slope_weight_2p0": (
        Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20_slope2.yaml"),
        Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20_slope2/pilot/best_model.pt"),
    ),
}


def summarize_split(result: dict[str, Any], probability_threshold: float) -> dict[str, Any]:
    sample_type = result["sample_type"]
    hard = sample_type == 1
    signal = sample_type == 2
    theta_true, _ = slopes_to_angles(result["slope_truth"])
    theta_pred, _ = slopes_to_angles(result["slope_prediction"])
    probability_pass = result["probabilities"] >= probability_threshold
    zero_hard_threshold = float(np.nextafter(result["probabilities"][hard].max(), 1.0))
    zero_hard_pass = result["probabilities"] >= zero_hard_threshold
    summary: dict[str, Any] = {
        "counts": {"hard": int(hard.sum()), "signal": int(signal.sum())},
        "hard_false_positive_at_p0p5": int((hard & probability_pass).sum()),
        "hard_false_positive_rate_at_p0p5": float(probability_pass[hard].mean()),
        "hard_probability_max": float(result["probabilities"][hard].max()),
        "zero_observed_hard_threshold": zero_hard_threshold,
        "signal_recall_at_zero_observed_hard_true_theta_ge_15": float(zero_hard_pass[signal & (theta_true >= 15)].mean()),
        "signal_recall_at_zero_observed_hard_true_theta_ge_20": float(zero_hard_pass[signal & (theta_true >= 20)].mean()),
        "signal_recall_p0p5_true_theta_ge_15": float(probability_pass[signal & (theta_true >= 15)].mean()),
        "signal_recall_p0p5_true_theta_ge_20": float(probability_pass[signal & (theta_true >= 20)].mean()),
        "selection_counts": {},
        "regression_by_true_theta": {},
    }
    for threshold in (15.0, 18.0, 20.0):
        pred_angle_pass = theta_pred >= threshold
        summary["selection_counts"][f"theta_pred_ge_{int(threshold)}"] = {
            "hard": int((hard & pred_angle_pass).sum()),
            "signal": int((signal & pred_angle_pass).sum()),
            "hard_and_p0p5": int((hard & pred_angle_pass & probability_pass).sum()),
            "signal_and_p0p5": int((signal & pred_angle_pass & probability_pass).sum()),
            "true_theta_ge_20_signal_efficiency": float(
                pred_angle_pass[signal & (theta_true >= 20.0)].mean()
            ),
            "true_theta_ge_20_signal_efficiency_and_p0p5": float(
                (pred_angle_pass & probability_pass)[signal & (theta_true >= 20.0)].mean()
            ),
        }
    for name, low, high in (("5_to_15", 5.0, 15.0), ("15_to_20", 15.0, 20.0), ("ge_20", 20.0, np.inf)):
        mask = signal & (theta_true >= low) & (theta_true < high)
        residual = theta_pred[mask] - theta_true[mask]
        relative = np.abs(residual) / theta_true[mask]
        summary["regression_by_true_theta"][name] = {
            "count": int(mask.sum()),
            "bias_mrad": float(residual.mean()) if mask.any() else None,
            "mae_mrad": float(np.abs(residual).mean()) if mask.any() else None,
            "mean_absolute_relative_error": float(relative.mean()) if mask.any() else None,
            "fraction_within_20_percent": float((relative <= 0.20).mean()) if mask.any() else None,
        }
    return summary


def main() -> None:
    report: dict[str, Any] = {}
    for name, (config_path, checkpoint_path) in PILOTS.items():
        config = load_config(config_path)
        device = torch.device("cpu")
        model = build_model(config).to(device)
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
        model.load_state_dict(checkpoint["model_state"])
        criterion = MultiTaskLoss(**config["loss"])
        data_config = CNNDataConfig(**config["data"]["dataset"], seed=int(config["seed"]))
        threshold = float(config["evaluation"]["threshold"])
        report[name] = {}
        for split_index, split in enumerate(("validation", "test"), start=1):
            loader = make_loader(
                Path(config["data"]["hdf5_path"]).resolve(), split, data_config,
                batch_size=int(config["training"]["batch_size"]), workers=0,
                seed=int(config["seed"]) + split_index,
                subset_size=config["subsets"]["pilot"][split], training=False, shuffle=False,
                include_sample_types=tuple(config["task"]["include_sample_types"]),
                positive_sample_types=tuple(config["task"]["positive_sample_types"]),
                signal_theta_range_mrad=None,
            )
            result = predict(
                model, loader, criterion, device,
                positive_sample_types=tuple(config["task"]["positive_sample_types"]),
                presence_positive_theta_min_mrad=float(config["task"]["presence_positive_theta_min_mrad"]),
            )
            report[name][split] = summarize_split(result, threshold)
    output = Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_pilot_comparison.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
