#!/usr/bin/env python3
"""Evaluate angle>15 classification stability on complete validation/test splits."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.metrics import average_precision_score, confusion_matrix, precision_recall_fscore_support, roc_auc_score

from training.cnn_dataset import CNNDataConfig
from training.cnn21d.losses import MultiTaskLoss
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, load_config, make_loader, predict


RUNS = {
    "seed_20260824": (
        Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20_slope2.yaml"),
        Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20_slope2/pilot/best_model.pt"),
    ),
    "seed_20260825": (
        Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20_slope2_seed20260825.yaml"),
        Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20_slope2_seed20260825/pilot/best_model.pt"),
    ),
    "seed_20260826": (
        Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20_slope2_seed20260826.yaml"),
        Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20_slope2_seed20260826/pilot/best_model.pt"),
    ),
}
OUTPUT_DIR = Path("runs/cnn21d_signal_t5_center40/angle15_seed_stability")


def classification_summary(labels: np.ndarray, probabilities: np.ndarray, threshold: float = 0.5) -> dict[str, Any]:
    prediction = probabilities >= threshold
    precision, recall, f1, _ = precision_recall_fscore_support(
        labels, prediction, average="binary", zero_division=0
    )
    return {
        "threshold": threshold,
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "auroc": float(roc_auc_score(labels, probabilities)),
        "auprc": float(average_precision_score(labels, probabilities)),
        "confusion_matrix": confusion_matrix(labels, prediction, labels=[0, 1]).tolist(),
    }


def summarize(result: dict[str, np.ndarray]) -> dict[str, Any]:
    sample_type = result["sample_type"]
    hard = sample_type == 1
    signal = sample_type == 2
    theta_true, _ = slopes_to_angles(result["slope_truth"])
    theta_pred, _ = slopes_to_angles(result["slope_prediction"])
    high_signal = signal & (theta_true > 15.0)
    labeled = hard | high_signal
    labels = high_signal[labeled].astype(np.int64)
    probabilities = result["probabilities"][labeled]
    selected = result["probabilities"] >= 0.5
    residual = theta_pred[signal] - theta_true[signal]
    output = {
        "counts": {
            "hard": int(hard.sum()),
            "signal_theta_gt_15": int(high_signal.sum()),
            "signal_theta_le_15_regression_only": int((signal & (theta_true <= 15.0)).sum()),
        },
        "classification": classification_summary(labels, probabilities),
        "hard_false_positives": int((hard & selected).sum()),
        "hard_false_positive_rate": float(selected[hard].mean()),
        "regression_all_signal": {
            "count": int(signal.sum()),
            "theta_mae_mrad": float(np.abs(residual).mean()),
            "theta_bias_mrad": float(residual.mean()),
            "theta_resolution_mrad": float(residual.std()),
        },
        "signal_recall_by_true_theta": {},
    }
    for name, low, high in (("15_to_20", 15.0, 20.0), ("20_to_30", 20.0, 30.0), ("ge_30", 30.0, np.inf)):
        mask = signal & (theta_true >= low) & (theta_true < high)
        output["signal_recall_by_true_theta"][name] = {
            "count": int(mask.sum()),
            "recall": float(selected[mask].mean()) if mask.any() else None,
            "mean_probability": float(result["probabilities"][mask].mean()) if mask.any() else None,
        }
    return output


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = {"task": {
        "bce_negative": "hard only",
        "bce_positive": "signal theta > 15 mrad",
        "bce_ignored": "signal theta <= 15 mrad",
        "regression": "all signal",
        "poisson": "excluded",
    }, "runs": {}}
    raw_by_split: dict[str, list[dict[str, np.ndarray]]] = {"validation": [], "test": []}

    for run_name, (config_path, checkpoint_path) in RUNS.items():
        config = load_config(config_path)
        model = build_model(config).cpu()
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
        model.load_state_dict(checkpoint["model_state"])
        criterion = MultiTaskLoss(**config["loss"])
        data_config = CNNDataConfig(**config["data"]["dataset"], seed=int(config["seed"]))
        report["runs"][run_name] = {"best_epoch": int(checkpoint["epoch"])}
        for split_index, split in enumerate(("validation", "test"), start=1):
            loader = make_loader(
                Path(config["data"]["hdf5_path"]).resolve(), split, data_config,
                batch_size=int(config["training"]["batch_size"]), workers=0,
                seed=int(config["seed"]) + split_index, subset_size=None,
                training=False, shuffle=False, include_sample_types=(1, 2),
                positive_sample_types=(2,), signal_theta_range_mrad=None,
            )
            result = predict(
                model, loader, criterion, torch.device("cpu"),
                positive_sample_types=(2,), presence_positive_theta_min_mrad=15.0,
            )
            raw_by_split[split].append(result)
            report["runs"][run_name][split] = summarize(result)

    report["ensemble"] = {}
    for split, results in raw_by_split.items():
        ensemble = dict(results[0])
        ensemble["probabilities"] = np.mean([result["probabilities"] for result in results], axis=0)
        ensemble["slope_prediction"] = np.mean([result["slope_prediction"] for result in results], axis=0)
        report["ensemble"][split] = summarize(ensemble)

    for split in ("validation", "test"):
        metrics = [report["runs"][name][split]["classification"] for name in RUNS]
        report.setdefault("seed_spread", {})[split] = {
            key: {"mean": float(np.mean([m[key] for m in metrics])), "std": float(np.std([m[key] for m in metrics]))}
            for key in ("precision", "recall", "f1", "auroc", "auprc")
        }

    (OUTPUT_DIR / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    names = list(RUNS) + ["ensemble"]
    x = np.arange(len(names))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, split in zip(axes, ("validation", "test"), strict=True):
        values = [report["runs"][name][split]["classification"]["f1"] for name in RUNS]
        values.append(report["ensemble"][split]["classification"]["f1"])
        ax.bar(x, values, color=["#4c72b0", "#55a868", "#c44e52", "#8172b3"])
        ax.set_xticks(x, ["seed 824", "seed 825", "seed 826", "ensemble"], rotation=20)
        ax.set_ylim(0.9, 1.005)
        ax.set_title(split)
        ax.grid(axis="y", alpha=0.25)
        for i, value in enumerate(values):
            ax.text(i, value + 0.002, f"{value:.3f}", ha="center", fontsize=9)
    axes[0].set_ylabel("F1 a soglia 0.5")
    fig.suptitle("Stabilità della classificazione θ > 15 mrad")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "f1_seed_stability.png", dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
