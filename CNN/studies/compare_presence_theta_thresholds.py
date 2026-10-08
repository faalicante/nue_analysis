#!/usr/bin/env python3
"""Compare presence checkpoints trained with different signal theta thresholds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from training.cnn21d.metrics import classification_metrics
from visualization.plot_pilot_diagnostics import collect_split
from training.train_cnn21d import build_model, load_config


def summarize_model(config_path: Path, checkpoint_path: Path) -> dict[str, object]:
    config = load_config(config_path)
    model = build_model(config).cpu()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    records = collect_split("test", config, model, torch.device("cpu"), mode="intermediate")
    signal = [row for row in records if row["sample_type"] == 2]
    hard = [row for row in records if row["sample_type"] == 1]
    bins: dict[str, object] = {}
    for name, low, high in (
        ("5_to_10", 5.0, 10.0),
        ("10_to_15", 10.0, 15.0),
        ("15_to_20", 15.0, 20.0),
        ("20_to_50", 20.0, 50.0),
        ("ge_50", 50.0, np.inf),
    ):
        selected = [row for row in signal if low <= row["theta_true_mrad"] < high]
        scores = np.asarray([row["probability"] for row in selected], dtype=np.float64)
        errors = np.asarray(
            [row["theta_pred_mrad"] - row["theta_true_mrad"] for row in selected],
            dtype=np.float64,
        )
        bins[name] = {
            "count": len(selected),
            "fraction_score_ge_0p5": float(np.mean(scores >= 0.5)) if scores.size else None,
            "score_median": float(np.median(scores)) if scores.size else None,
            "score_q10": float(np.quantile(scores, 0.1)) if scores.size else None,
            "theta_mae_mrad": float(np.mean(np.abs(errors))) if errors.size else None,
        }

    hard_scores = np.asarray([row["probability"] for row in hard], dtype=np.float64)
    common_tasks: dict[str, object] = {}
    for cutoff in (10.0, 15.0):
        eligible_signal = [row for row in signal if row["theta_true_mrad"] > cutoff]
        labels = np.concatenate((np.zeros(len(hard)), np.ones(len(eligible_signal))))
        probabilities = np.concatenate((
            hard_scores,
            np.asarray([row["probability"] for row in eligible_signal], dtype=np.float64),
        ))
        common_tasks[f"signal_gt_{int(cutoff)}"] = classification_metrics(
            labels, probabilities, threshold=0.5
        )
    return {
        "checkpoint_epoch": int(checkpoint["epoch"]),
        "trained_presence_minimum_mrad": config["task"]["presence_positive_theta_min_mrad"],
        "angle_bins": bins,
        "hard": {
            "count": len(hard),
            "fraction_score_ge_0p5": float(np.mean(hard_scores >= 0.5)),
            "score_median": float(np.median(hard_scores)),
            "score_q90": float(np.quantile(hard_scores, 0.9)),
            "score_q99": float(np.quantile(hard_scores, 0.99)),
            "score_max": float(hard_scores.max()),
        },
        "common_evaluation_tasks": common_tasks,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--model", action="append", nargs=3, metavar=("NAME", "CONFIG", "CHECKPOINT"), required=True
    )
    args = parser.parse_args()
    result = {
        name: summarize_model(Path(config), Path(checkpoint))
        for name, config, checkpoint in args.model
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
