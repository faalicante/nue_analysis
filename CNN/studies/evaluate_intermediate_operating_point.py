#!/usr/bin/env python3
"""Choose a presence threshold on validation and apply it once to test."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import precision_recall_curve

from training.cnn_dataset import CNNDataConfig
from training.cnn21d.losses import MultiTaskLoss
from training.cnn21d.metrics import classification_metrics, slopes_to_angles
from training.train_cnn21d import build_model, load_config, make_loader, predict


def infer(
    config_path: Path, checkpoint_path: Path, split: str, split_index: int
) -> dict[str, np.ndarray]:
    config = load_config(config_path)
    model = build_model(config).cpu()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    criterion = MultiTaskLoss(**config["loss"])
    data_config = CNNDataConfig(**config["data"]["dataset"], seed=int(config["seed"]))
    loader = make_loader(
        Path(config["data"]["hdf5_path"]).resolve(), split, data_config,
        batch_size=int(config["training"]["batch_size"]), workers=0,
        seed=int(config["seed"]) + split_index, subset_size=None,
        training=False, shuffle=False, include_sample_types=(1, 2),
        positive_sample_types=(2,), signal_theta_range_mrad=None,
    )
    return predict(
        model, loader, criterion, torch.device("cpu"), positive_sample_types=(2,),
        presence_positive_theta_min_mrad=15.0,
    )


def labeled_arrays(result: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    mask = result["classification_mask"].astype(bool)
    return result["labels"][mask].astype(np.int64), result["probabilities"][mask]


def per_angle_recall(result: dict[str, np.ndarray], threshold: float) -> dict[str, dict[str, float | int]]:
    signal = result["sample_type"] == 2
    theta, _ = slopes_to_angles(result["slope_truth"])
    selected = result["probabilities"] >= threshold
    output = {}
    for name, low, high in (("15_to_20", 15.0, 20.0), ("20_to_30", 20.0, 30.0), ("ge_30", 30.0, np.inf)):
        mask = signal & (theta > low) & (theta < high)
        output[name] = {"count": int(mask.sum()), "recall": float(selected[mask].mean())}
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path,
        default=Path("training/configs/cnn21d_signal_gt15_hard_only_intermediate4000.yaml"),
    )
    parser.add_argument(
        "--checkpoint", type=Path,
        default=Path(
            "runs/cnn21d_signal_t5_center40/signal_gt15_hard_only_intermediate4000/"
            "intermediate/best_model.pt"
        ),
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    output = args.output or args.checkpoint.parent / "operating_point.json"
    validation = infer(args.config, args.checkpoint, "validation", 1)
    test = infer(args.config, args.checkpoint, "test", 2)
    y_val, p_val = labeled_arrays(validation)
    precision, recall, thresholds = precision_recall_curve(y_val, p_val)
    f1 = 2.0 * precision[:-1] * recall[:-1] / np.maximum(precision[:-1] + recall[:-1], 1e-12)
    best_index = int(np.argmax(f1))
    threshold = float(thresholds[best_index])
    operating_points = {"max_validation_f1": threshold}
    for target_recall in (0.99, 0.995):
        eligible = np.flatnonzero(recall[:-1] >= target_recall)
        operating_points[f"highest_threshold_with_validation_recall_ge_{str(target_recall).replace('.', 'p')}"] = (
            float(thresholds[eligible[-1]])
        )
    report = {
        "selection_rule": "all thresholds selected on validation only",
        "operating_points": {
            name: {
                "threshold": selected_threshold,
                "validation": classification_metrics(y_val, p_val, selected_threshold),
                "test": classification_metrics(*labeled_arrays(test), selected_threshold),
                "test_recall_by_true_theta": per_angle_recall(test, selected_threshold),
            }
            for name, selected_threshold in operating_points.items()
        },
        "fixed_threshold_0p5": {
            "validation": classification_metrics(y_val, p_val, 0.5),
            "test": classification_metrics(*labeled_arrays(test), 0.5),
        },
    }
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
