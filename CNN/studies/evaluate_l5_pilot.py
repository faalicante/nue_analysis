#!/usr/bin/env python3
"""Compare the low-angle pilot with the current full model on held-out data."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import torch

from training.cnn_dataset import CNNDataConfig, normalize_per_slice
from training.cnn21d.losses import MultiTaskLoss
from training.train_cnn21d import (
    build_model, combine_head_predictions, load_config, make_loader, predict,
    result_metrics,
)


ANGLE_BINS = (
    (0.0, 5.0), (5.0, 10.0), (10.0, 20.0), (20.0, 50.0),
    (50.0, 100.0), (100.0, 200.0),
)


def theta_mrad(slopes: np.ndarray) -> np.ndarray:
    return 1000.0 * np.arctan(np.linalg.norm(slopes, axis=1) / 27.0)


def load_checkpoint_model(path: Path, fallback_config: dict[str, Any]) -> torch.nn.Module:
    checkpoint = torch.load(path, map_location="cpu", weights_only=False)
    config = checkpoint.get("config", fallback_config)
    model = build_model(config)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    return model


def signal_angle_groups(result: dict[str, Any]) -> dict[str, Any]:
    truth_theta = theta_mrad(result["slope_truth"])
    predicted_theta = theta_mrad(result["slope_prediction"])
    signal = result["sample_type"] == 2
    groups: dict[str, Any] = {}
    for low, high in ANGLE_BINS:
        mask = signal & (truth_theta > low) & (truth_theta <= high)
        error = predicted_theta[mask] - truth_theta[mask]
        groups[f"{low:g}-{high:g}"] = {
            "count": int(mask.sum()),
            "mean_score": float(result["probabilities"][mask].mean()) if mask.any() else None,
            "score_gt_0p5": float(np.mean(result["probabilities"][mask] > 0.5)) if mask.any() else None,
            "score_gt_0p9": float(np.mean(result["probabilities"][mask] > 0.9)) if mask.any() else None,
            "theta_pred_mean_mrad": float(predicted_theta[mask].mean()) if mask.any() else None,
            "theta_bias_mrad": float(error.mean()) if mask.any() else None,
            "theta_mae_mrad": float(np.abs(error).mean()) if mask.any() else None,
        }
    return groups


def negative_groups(result: dict[str, Any]) -> dict[str, Any]:
    groups: dict[str, Any] = {}
    for code, name in ((0, "poisson"), (1, "hard")):
        mask = result["sample_type"] == code
        groups[name] = {
            "count": int(mask.sum()),
            "mean_score": float(result["probabilities"][mask].mean()),
            "score_gt_0p5": float(np.mean(result["probabilities"][mask] > 0.5)),
            "score_gt_0p9": float(np.mean(result["probabilities"][mask] > 0.9)),
            "score_gt_0p98": float(np.mean(result["probabilities"][mask] > 0.98)),
        }
    return groups


def predict_ranked_data_crops(
    manifest_path: Path,
    models: dict[str, tuple[torch.nn.Module, torch.nn.Module]],
) -> list[dict[str, Any]]:
    candidates = json.loads(manifest_path.read_text(encoding="utf-8"))["candidates"]
    rows: list[dict[str, Any]] = []
    for candidate in candidates:
        with h5py.File(candidate["volume_path"], "r") as source:
            event_index = int(candidate["event_index"])
            x0 = int(candidate["x_start"])
            y0 = int(candidate["y_start"])
            size = int(candidate["crop_size"])
            raw = np.asarray(
                source["volumes_raw"][event_index, :, y0:y0 + size, x0:x0 + size]
            )
            mu = float(source["background_mu"][event_index])
        volume = torch.from_numpy(normalize_per_slice(raw, mu)[None, None]).float()
        row: dict[str, Any] = {
            "rank": int(candidate["rank"]),
            "cell_id": int(candidate["cell_id"]),
            "old_scan_score": float(candidate["presence_score"]),
            "old_scan_theta_mrad": float(candidate["theta_pred_mrad"]),
        }
        for name, (classification_model, regression_model) in models.items():
            with torch.no_grad():
                classification_output = classification_model(volume)
                regression_output = regression_model(volume)
            slope = regression_output["slope_xy"].numpy()
            row[f"{name}_score"] = float(torch.sigmoid(classification_output["presence_logit"])[0])
            row[f"{name}_theta_mrad"] = float(theta_mrad(slope)[0])
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--new-checkpoint", type=Path)
    parser.add_argument("--old-checkpoint", type=Path)
    parser.add_argument("--control-checkpoint", type=Path)
    parser.add_argument("--poisson-checkpoint", type=Path)
    parser.add_argument("--dual-classification-checkpoint", type=Path)
    parser.add_argument("--dual-regression-checkpoint", type=Path)
    parser.add_argument(
        "--layer-drop-model",
        choices=("new_l5", "old_full", "control", "with_poisson", "dual_balanced"),
        help="Model on which to run the full-validation layer-drop diagnostic.",
    )
    parser.add_argument("--data-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    path = Path(config["data"]["hdf5_path"]).resolve()
    data_cfg = CNNDataConfig(**config["data"]["dataset"], seed=int(config["seed"]))
    criterion = MultiTaskLoss(**config["loss"])
    checkpoints = {
        "new_l5": args.new_checkpoint,
        "old_full": args.old_checkpoint,
        "control": args.control_checkpoint,
        "with_poisson": args.poisson_checkpoint,
    }
    single_models = {
        name: load_checkpoint_model(checkpoint, config)
        for name, checkpoint in checkpoints.items()
        if checkpoint is not None
    }
    models = {name: (model, model) for name, model in single_models.items()}
    if bool(args.dual_classification_checkpoint) != bool(args.dual_regression_checkpoint):
        parser.error("dual classification and regression checkpoints must be provided together")
    if args.dual_classification_checkpoint is not None:
        models["dual_balanced"] = (
            load_checkpoint_model(args.dual_classification_checkpoint, config),
            load_checkpoint_model(args.dual_regression_checkpoint, config),
        )
    if not models:
        parser.error("provide at least one checkpoint")
    if args.layer_drop_model is not None and args.layer_drop_model not in models:
        parser.error(f"--layer-drop-model {args.layer_drop_model!r} has no checkpoint")

    summaries: dict[str, Any] = {}
    cached: dict[tuple[str, str], dict[str, Any]] = {}
    for model_name, (classification_model, regression_model) in models.items():
        summaries[model_name] = {}
        for split_index, split in enumerate(("validation", "test"), start=1):
            loader = make_loader(
                path, split, data_cfg,
                batch_size=32, workers=0, seed=int(config["seed"]) + split_index,
                subset_size=None, training=False, shuffle=False,
                include_sample_types=(0, 1, 2), positive_sample_types=(2,),
                signal_theta_range_mrad=None,
            )
            classification_result = predict(
                classification_model, loader, criterion, torch.device("cpu"),
                positive_sample_types=(2,),
                presence_positive_theta_min_mrad=10.0,
                presence_low_angle_policy="negative",
            )
            if regression_model is classification_model:
                result = classification_result
            else:
                regression_result = predict(
                    regression_model, loader, criterion, torch.device("cpu"),
                    positive_sample_types=(2,),
                    presence_positive_theta_min_mrad=10.0,
                    presence_low_angle_policy="negative",
                )
                result = combine_head_predictions(classification_result, regression_result)
            cached[(model_name, split)] = result
            summaries[model_name][split] = {
                "all": result_metrics(result, 0.5),
                "signal_angle_bins": signal_angle_groups(result),
                "background": negative_groups(result),
            }

    layer_drop: dict[str, Any] = {}
    if args.layer_drop_model is not None:
        baseline = cached[(args.layer_drop_model, "validation")]
        classification_model, regression_model = models[args.layer_drop_model]
        for transform in (
            "drop_one_random", "drop_two_isolated", "drop_two_consecutive", "drop_peak"
        ):
            loader = make_loader(
                path, "validation", data_cfg,
                batch_size=32, workers=0, seed=int(config["seed"]) + 1,
                subset_size=None, training=False, shuffle=False,
                include_sample_types=(0, 1, 2), positive_sample_types=(2,),
                signal_theta_range_mrad=None,
            )
            classification_result = predict(
                classification_model, loader, criterion, torch.device("cpu"),
                z_transform=transform, z_seed=int(config["seed"]),
                positive_sample_types=(2,), presence_positive_theta_min_mrad=10.0,
                presence_low_angle_policy="negative",
            )
            if regression_model is classification_model:
                result = classification_result
            else:
                regression_result = predict(
                    regression_model, loader, criterion, torch.device("cpu"),
                    z_transform=transform, z_seed=int(config["seed"]),
                    positive_sample_types=(2,), presence_positive_theta_min_mrad=10.0,
                    presence_low_angle_policy="negative",
                )
                result = combine_head_predictions(classification_result, regression_result)
            layer_drop[transform] = {
                "all": result_metrics(result, 0.5),
                "mean_probability_abs_change": float(np.abs(
                    baseline["probabilities"] - result["probabilities"]
                ).mean()),
                "signal_angle_bins": signal_angle_groups(result),
            }

    data_rows = predict_ranked_data_crops(args.data_manifest, models)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = {
        "dataset": str(path),
        "models": summaries,
        "layer_drop_validation": {
            "model": args.layer_drop_model,
            "results": layer_drop,
        },
        "data_scan_candidates": data_rows,
    }
    (args.output_dir / "evaluation.json").write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    with (args.output_dir / "data_candidate_rescoring.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data_rows[0]))
        writer.writeheader()
        writer.writerows(data_rows)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
