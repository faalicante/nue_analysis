#!/usr/bin/env python3
"""Report pilot performance in signal theta bins for the preselected t5 sample."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

from visualization.plot_pilot_diagnostics import collect_split
from training.train_cnn21d import build_model, load_config


def summarize(config_path: Path, checkpoint_path: Path) -> dict[str, object]:
    config = load_config(config_path)
    model = build_model(config).cpu()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    records = collect_split("validation", config, model, torch.device("cpu")) + collect_split("test", config, model, torch.device("cpu"))
    result: dict[str, object] = {}
    for split in ("validation", "test"):
        signal = [row for row in records if row["split"] == split and row["sample_type"] == 2]
        split_result: dict[str, object] = {}
        for name, low, high in (("5_to_10", 5.0, 10.0), ("10_to_50", 10.0, 50.0), ("ge_50", 50.0, np.inf)):
            selected = [row for row in signal if low <= row["theta_true_mrad"] < high]
            errors = np.asarray([row["theta_pred_mrad"] - row["theta_true_mrad"] for row in selected])
            split_result[name] = {
                "count": len(selected),
                "classification_recall": float(np.mean([row["predicted_label"] for row in selected])) if selected else None,
                "theta_bias_mrad": float(errors.mean()) if len(errors) else None,
                "theta_mae_mrad": float(np.abs(errors).mean()) if len(errors) else None,
                "theta_resolution_mrad": float(errors.std(ddof=0)) if len(errors) else None,
            }
        result[split] = split_result
    return result


def main() -> None:
    output: dict[str, object] = {}
    for crop in (20, 32):
        base = Path(f"runs/cnn21d_signal_t5_center40/all_t5/crop{crop}/pilot")
        output[f"crop{crop}"] = summarize(
            Path(f"training/configs/cnn21d_signal_t5_center40_all_crop{crop}.yaml"), base / "best_model.pt"
        )
    path = Path("runs/cnn21d_signal_t5_center40/all_t5/signal_angle_bin_metrics.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
