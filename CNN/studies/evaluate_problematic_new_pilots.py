#!/usr/bin/env python3
"""Evaluate the previously problematic signal events with the new centering pilots."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import torch

from training.cnn_dataset import extract_xy_crop, normalize_per_slice
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, load_config


DEFAULT_EVENT_IDS = (5461, 3736, 938, 9255, 5022, 6008)


def load_pilot(config_path: Path, checkpoint_path: Path) -> tuple[dict[str, Any], torch.nn.Module]:
    config = load_config(config_path)
    model = build_model(config).cpu()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    return config, model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hdf5", type=Path, default=Path("cnn_dataset_signal_t5_center40.h5"))
    parser.add_argument("--output", type=Path, default=Path("runs/cnn21d_signal_t5_center40/problematic_predictions.json"))
    parser.add_argument("--event-id", type=int, action="append", dest="event_ids")
    args = parser.parse_args()
    event_ids = tuple(args.event_ids or DEFAULT_EVENT_IDS)
    pilots = {
        "crop20_all_t5": (
            Path("training/configs/cnn21d_signal_t5_center40_all_crop20.yaml"),
            Path("runs/cnn21d_signal_t5_center40/all_t5/crop20/pilot/best_model.pt"),
        ),
        "crop32_all_t5": (
            Path("training/configs/cnn21d_signal_t5_center40_all_crop32.yaml"),
            Path("runs/cnn21d_signal_t5_center40/all_t5/crop32/pilot/best_model.pt"),
        ),
    }
    loaded = {name: load_pilot(*paths) for name, paths in pilots.items()}
    rows: list[dict[str, Any]] = []
    with h5py.File(args.hdf5, "r") as hdf5:
        ids = np.asarray(hdf5["signal_event_id"][:], dtype=np.int64)
        sample_types = np.asarray(hdf5["sample_type"][:], dtype=np.int64)
        for event_id in event_ids:
            matches = np.flatnonzero((ids == event_id) & (sample_types == 2))
            if len(matches) != 1:
                rows.append({"event_id": event_id, "status": f"found_{len(matches)}"})
                continue
            index = int(matches[0])
            raw = np.asarray(hdf5["regions_raw"][index], dtype=np.float32)
            background_values = np.asarray(hdf5["background_mu"][index], dtype=np.float32)
            mu = float(background_values.reshape(-1)[0])
            truth = np.asarray(hdf5["slope_xy"][index], dtype=np.float32)
            theta_true = float(slopes_to_angles(truth[None])[0][0])
            row: dict[str, Any] = {
                "event_id": event_id,
                "hdf5_index": index,
                "split": ("train", "validation", "test")[int(hdf5["split"][index])],
                "slope_true": truth.tolist(),
                "theta_true_mrad": theta_true,
                "pilots": {},
            }
            for name, (config, model) in loaded.items():
                crop_size = int(config["data"]["dataset"]["crop_size"])
                normalized = normalize_per_slice(extract_xy_crop(raw, crop_size=crop_size), mu)
                volume = torch.from_numpy(normalized[None, None])
                with torch.no_grad():
                    output = model(volume)
                probability = float(torch.sigmoid(output["presence_logit"])[0])
                prediction = output["slope_xy"][0].numpy()
                theta_pred = float(slopes_to_angles(prediction[None])[0][0])
                row["pilots"][name] = {
                    "probability_signal": probability,
                    "predicted_signal": probability >= float(config["evaluation"]["threshold"]),
                    "slope_pred": prediction.tolist(),
                    "theta_pred_mrad": theta_pred,
                    "theta_abs_error_mrad": abs(theta_pred - theta_true),
                }
            rows.append(row)
    result = {"hdf5": str(args.hdf5.resolve()), "events": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
