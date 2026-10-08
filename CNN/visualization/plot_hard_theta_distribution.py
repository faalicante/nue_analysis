#!/usr/bin/env python3
"""Plot the predicted-angle distribution for the unbiased hard control sample."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch

from training.cnn_dataset import CNNDataConfig
from training.cnn21d.losses import MultiTaskLoss
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, load_config, make_loader, predict


CONFIG = Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20_slope2.yaml")
CHECKPOINT = Path(
    "runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20_slope2/pilot/best_model.pt"
)
OUTPUT_DIR = Path("runs/cnn21d_signal_t5_center40/hard_theta_pred")


def main() -> None:
    config = load_config(CONFIG)
    model = build_model(config).cpu()
    checkpoint = torch.load(CHECKPOINT, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    criterion = MultiTaskLoss(**config["loss"])
    data_config = CNNDataConfig(**config["data"]["dataset"], seed=int(config["seed"]))

    results = []
    split_labels = []
    for split_index, split in enumerate(("validation", "test"), start=1):
        loader = make_loader(
            Path(config["data"]["hdf5_path"]).resolve(),
            split,
            data_config,
            batch_size=int(config["training"]["batch_size"]),
            workers=0,
            seed=int(config["seed"]) + split_index,
            subset_size=None,
            training=False,
            shuffle=False,
            include_sample_types=(1,),
            positive_sample_types=(2,),
            signal_theta_range_mrad=None,
        )
        result = predict(
            model,
            loader,
            criterion,
            torch.device("cpu"),
            positive_sample_types=(2,),
            presence_positive_theta_min_mrad=15.0,
        )
        results.append(result)
        split_labels.extend([split] * len(result["probabilities"]))

    slope_prediction = np.concatenate([result["slope_prediction"] for result in results])
    probability = np.concatenate([result["probabilities"] for result in results])
    theta_pred, _ = slopes_to_angles(slope_prediction)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUTPUT_DIR / "hard_theta_pred.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("split", "theta_pred_mrad", "presence_probability"))
        writer.writerows(zip(split_labels, theta_pred, probability, strict=True))

    upper = max(30.0, float(np.ceil(theta_pred.max() / 2.0) * 2.0))
    bins = np.arange(0.0, upper + 0.5, 0.5)
    counts, edges = np.histogram(theta_pred, bins=bins)
    summary = {
        "sample": "hard validation+test, no jitter",
        "count": int(theta_pred.size),
        "mean_mrad": float(theta_pred.mean()),
        "median_mrad": float(np.median(theta_pred)),
        "std_mrad": float(theta_pred.std()),
        "quantiles_mrad": {
            str(q): float(np.quantile(theta_pred, q))
            for q in (0.5, 0.9, 0.95, 0.99, 0.995, 0.999)
        },
        "counts_above_mrad": {
            str(threshold): int((theta_pred >= threshold).sum())
            for threshold in (5, 10, 15, 20, 25, 30)
        },
        "counts_above_mrad_and_presence_ge_0p5": {
            str(threshold): int(((theta_pred >= threshold) & (probability >= 0.5)).sum())
            for threshold in (5, 10, 15, 20, 25, 30)
        },
        "histogram": {
            "bin_edges_mrad": edges.tolist(),
            "counts": counts.tolist(),
        },
    }
    (OUTPUT_DIR / "hard_theta_pred_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(theta_pred, bins=bins, histtype="stepfilled", alpha=0.72, color="#3366cc")
    ax.axvline(15, color="#e68600", linestyle="--", linewidth=1.8, label="15 mrad")
    ax.axvline(20, color="#c62828", linestyle="--", linewidth=1.8, label="20 mrad")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\theta_{pred}$ [mrad]")
    ax.set_ylabel("Eventi hard / 0.5 mrad")
    ax.set_title(r"Distribuzione di $\theta_{pred}$ — hard validation + test")
    ax.grid(alpha=0.22, which="both")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "hard_theta_pred_distribution.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
