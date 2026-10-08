#!/usr/bin/env python3
"""Estimate candidate yield by folding pilot response with cumulative angle efficiencies."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

from training.cnn_dataset import CNNDataConfig
from training.cnn21d.losses import MultiTaskLoss
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, load_config, make_loader, predict


SIGNAL_CUMULATIVE_PERCENT = {
    0: 85.51, 1: 85.14, 2: 83.69, 3: 81.0, 4: 77.79, 5: 74.28,
    6: 70.99, 7: 67.23, 8: 63.64, 9: 59.96, 10: 56.73,
    11: 53.37, 12: 50.02, 13: 47.23, 14: 44.58, 15: 41.86,
    16: 39.21, 17: 36.91, 18: 34.8, 19: 32.96, 20: 31.06,
}
MUON_CUMULATIVE_PERCENT = {
    0: 100.0, 1: 88.12, 2: 62.51, 3: 39.17, 4: 25.15, 5: 16.07,
    6: 8.7, 7: 4.31, 8: 2.6, 9: 1.89, 10: 1.36, 11: 0.92,
    12: 0.61, 13: 0.42, 14: 0.3, 15: 0.23, 16: 0.17, 17: 0.12,
    18: 0.08, 19: 0.04, 20: 0.01,
}


def main() -> None:
    config = load_config(Path("training/configs/cnn21d_signal_t5_center40_angle15_jitter5_crop20_slope2.yaml"))
    model = build_model(config).cpu()
    checkpoint = torch.load(
        "runs/cnn21d_signal_t5_center40/angle15_jitter5_crop20_slope2/pilot/best_model.pt",
        map_location="cpu", weights_only=False,
    )
    model.load_state_dict(checkpoint["model_state"])
    criterion = MultiTaskLoss(**config["loss"])
    data_config = CNNDataConfig(**config["data"]["dataset"], seed=int(config["seed"]))
    collected = []
    for split_index, split in enumerate(("validation", "test"), start=1):
        loader = make_loader(
            Path(config["data"]["hdf5_path"]).resolve(), split, data_config,
            batch_size=int(config["training"]["batch_size"]), workers=0,
            seed=int(config["seed"]) + split_index,
            subset_size=None, training=False, shuffle=False,
            include_sample_types=(1, 2), positive_sample_types=(2,), signal_theta_range_mrad=None,
        )
        collected.append(predict(
            model, loader, criterion, torch.device("cpu"), positive_sample_types=(2,),
            presence_positive_theta_min_mrad=15.0,
        ))
    sample_type = np.concatenate([item["sample_type"] for item in collected])
    probability = np.concatenate([item["probabilities"] for item in collected])
    slope_truth = np.concatenate([item["slope_truth"] for item in collected])
    slope_prediction = np.concatenate([item["slope_prediction"] for item in collected])
    theta_true, _ = slopes_to_angles(slope_truth)
    theta_pred, _ = slopes_to_angles(slope_prediction)
    signal = sample_type == 2
    hard = sample_type == 1

    def fold(probability_threshold: float) -> tuple[float, float]:
        selected_at_threshold = (probability >= probability_threshold) & (theta_pred >= 20.0)
        muon_candidates = 0.0
        signal_efficiency_percent = 0.0
        for low in range(5, 20):
            mask = signal & (theta_true >= low) & (theta_true < low + 1)
            if not mask.any():
                continue
            response = float(selected_at_threshold[mask].mean())
            muon_fraction = (MUON_CUMULATIVE_PERCENT[low] - MUON_CUMULATIVE_PERCENT[low + 1]) / 100.0
            signal_fraction = (SIGNAL_CUMULATIVE_PERCENT[low] - SIGNAL_CUMULATIVE_PERCENT[low + 1]) / 100.0
            muon_candidates += 1.0e6 * muon_fraction * response
            signal_efficiency_percent += 100.0 * signal_fraction * response
        mask = signal & (theta_true >= 20.0)
        response = float(selected_at_threshold[mask].mean())
        muon_candidates += 1.0e6 * MUON_CUMULATIVE_PERCENT[20] / 100.0 * response
        signal_efficiency_percent += SIGNAL_CUMULATIVE_PERCENT[20] * response
        return muon_candidates, signal_efficiency_percent

    selected = (probability >= 0.5) & (theta_pred >= 20.0)

    bins = []
    expected_muon_candidates = 0.0
    expected_signal_efficiency_percent = 0.0
    for low in range(5, 20):
        mask = signal & (theta_true >= low) & (theta_true < low + 1)
        response = float(selected[mask].mean()) if mask.any() else None
        muon_fraction = (MUON_CUMULATIVE_PERCENT[low] - MUON_CUMULATIVE_PERCENT[low + 1]) / 100.0
        signal_fraction = (SIGNAL_CUMULATIVE_PERCENT[low] - SIGNAL_CUMULATIVE_PERCENT[low + 1]) / 100.0
        contribution = None if response is None else 1.0e6 * muon_fraction * response
        if contribution is not None:
            expected_muon_candidates += contribution
            expected_signal_efficiency_percent += 100.0 * signal_fraction * response
        bins.append({"theta_bin_mrad": [low, low + 1], "pilot_signal_count": int(mask.sum()),
                     "selection_response": response, "expected_muon_candidates": contribution})
    mask = signal & (theta_true >= 20.0)
    response = float(selected[mask].mean())
    high_contribution = 1.0e6 * MUON_CUMULATIVE_PERCENT[20] / 100.0 * response
    expected_muon_candidates += high_contribution
    expected_signal_efficiency_percent += SIGNAL_CUMULATIVE_PERCENT[20] * response
    bins.append({"theta_bin_mrad": [20, None], "pilot_signal_count": int(mask.sum()),
                 "selection_response": response, "expected_muon_candidates": high_contribution})

    result = {
        "selection": "presence_probability >= 0.5 and theta_pred >= 20 mrad",
        "assumed_showering_muons_before_angle_selection": 1_000_000,
        "estimated_muon_candidates_from_true_theta_ge_5": expected_muon_candidates,
        "estimated_signal_efficiency_percent_of_original": expected_signal_efficiency_percent,
        "presence_threshold_scan_at_theta_pred_ge_20": [
            {
                "presence_threshold": threshold,
                "estimated_muon_candidates": fold(threshold)[0],
                "estimated_signal_efficiency_percent_of_original": fold(threshold)[1],
            }
            for threshold in (0.5, 0.9, 0.95, 0.975, 0.99, 0.995, 0.996, 0.997, 0.998, 0.999)
        ],
        "hard_low_angle_control": {
            "count": int(hard.sum()), "selected": int(selected[hard].sum()),
            "note": "Zero observed does not constrain a 1e-4 tail with this sample size.",
        },
        "bins": bins,
        "limitations": [
            "Response is estimated from the full validation+test signal samples; per-angle bins still have finite statistics.",
            "Bins with no pilot signal cannot be folded and are omitted.",
            "The estimate assumes identical CNN response for electron and muon showers at fixed true theta.",
            "Overlapping-window duplicates must be removed before counting manual candidates.",
        ],
    }
    output = Path("runs/cnn21d_signal_t5_center40/angle15_jitter5_candidate_yield.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
