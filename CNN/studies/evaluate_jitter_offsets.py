#!/usr/bin/env python3
"""Evaluate the angle15 crop20 pilot at fixed scan-grid miscenterings."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset

from training.cnn_dataset import CNNDataConfig, HDF5CNN21DDataset
from training.cnn21d.losses import MultiTaskLoss
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import _stratified_local_indices, build_model, load_config, predict, result_metrics


class FixedOffsetDataset(HDF5CNN21DDataset):
    def __init__(self, *args, offset_xy: tuple[int, int], **kwargs):
        self.fixed_offset_xy = offset_xy
        super().__init__(*args, **kwargs)

    def _augmentation_parameters(self, raw_index: int, is_shower_class: bool):
        del raw_index, is_shower_class
        return self.fixed_offset_xy[0], self.fixed_offset_xy[1], 0, False, False


def main() -> None:
    config_path = Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20_slope2.yaml")
    checkpoint_path = Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20_slope2/pilot/best_model.pt")
    config = load_config(config_path)
    model = build_model(config).cpu()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    criterion = MultiTaskLoss(**config["loss"])
    data_config = CNNDataConfig(**config["data"]["dataset"], seed=int(config["seed"]))
    report = {}
    for offset_x, offset_y in ((0, 0), (-5, 0), (5, 0), (0, -5), (0, 5), (-5, -5), (-5, 5), (5, -5), (5, 5)):
        dataset = FixedOffsetDataset(
            Path(config["data"]["hdf5_path"]), split="validation", training=False,
            config=data_config, offset_xy=(offset_x, offset_y),
        )
        indices = _stratified_local_indices(
            dataset, int(config["subsets"]["pilot"]["validation"]), int(config["seed"]) + 1,
            tuple(config["task"]["include_sample_types"]), tuple(config["task"]["positive_sample_types"]), None,
        )
        loader = DataLoader(Subset(dataset, indices), batch_size=int(config["training"]["batch_size"]), shuffle=False)
        result = predict(
            model, loader, criterion, torch.device("cpu"),
            positive_sample_types=tuple(config["task"]["positive_sample_types"]),
            presence_positive_theta_min_mrad=float(config["task"]["presence_positive_theta_min_mrad"]),
        )
        metrics = result_metrics(result, float(config["evaluation"]["threshold"]))
        theta_true, _ = slopes_to_angles(result["slope_truth"])
        theta_pred, _ = slopes_to_angles(result["slope_prediction"])
        high = (result["sample_type"] == 2) & (theta_true >= 20.0)
        residual = theta_pred[high] - theta_true[high]
        report[f"x{offset_x:+d}_y{offset_y:+d}"] = {
            "classification": {key: metrics[key] for key in ("precision", "recall", "f1", "auroc", "auprc", "confusion_matrix")},
            "theta_ge_20": {
                "count": int(high.sum()),
                "mae_mrad": float(np.abs(residual).mean()),
                "bias_mrad": float(residual.mean()),
                "fraction_within_20_percent": float((np.abs(residual) / theta_true[high] <= 0.20).mean()),
            },
        }
    output = Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_fixed_offset_validation.json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
