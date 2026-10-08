#!/usr/bin/env python3
"""Compare z-order sensitivity of scan candidates and matched signal controls."""

from __future__ import annotations

import glob
import json
from pathlib import Path

import h5py
import numpy as np
import torch

from training.cnn_dataset import normalize_per_slice
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, load_config


CONFIG = Path("training/configs/cnn21d_signal_gt10_hard_only_intermediate4000.yaml")
CHECKPOINT = Path(
    "runs/cnn21d_signal_t5_center40/"
    "signal_gt10_hard_only_intermediate4000/intermediate/best_model.pt"
)
BACKGROUND_VOLUMES = "/Users/fabioali/cernbox/CNN/background_scan_volumes/background_scan_*.h5"
BACKGROUND_PREDICTIONS = (
    "runs/cnn21d_signal_t5_center40/background_scan_full_t098/"
    "background_scan_*_t098.h5"
)
SIGNAL_VOLUMES = Path("/private/tmp/crop_bw_scan_splits/signal_scan_test.h5")
SIGNAL_PREDICTIONS = Path(
    "runs/cnn21d_signal_t5_center40/scan_volumes_test_gt10_t098.h5"
)


def collect_background() -> np.ndarray:
    sources = {path.name: path for path in map(Path, glob.glob(BACKGROUND_VOLUMES))}
    crops: list[np.ndarray] = []
    for prediction_path in sorted(map(Path, glob.glob(BACKGROUND_PREDICTIONS))):
        source_name = prediction_path.name.removesuffix("_t098.h5") + ".h5"
        with h5py.File(sources[source_name], "r") as source, h5py.File(prediction_path, "r") as prediction:
            group = prediction["candidates"]
            xs = prediction["grid_x_start_bin"][:]
            ys = prediction["grid_y_start_bin"][:]
            size = int(prediction.attrs["crop_size"])
            for i in range(len(group["event_id"])):
                event = int(group["event_index"][i])
                x0 = int(xs[int(group["grid_col"][i])])
                y0 = int(ys[int(group["grid_row"][i])])
                raw = source["volumes_raw"][event, :, y0:y0 + size, x0:x0 + size]
                crops.append(normalize_per_slice(raw, source["background_mu"][event]))
    return np.stack(crops)


def collect_signal() -> np.ndarray:
    crops: list[np.ndarray] = []
    with h5py.File(SIGNAL_VOLUMES, "r") as source, h5py.File(SIGNAL_PREDICTIONS, "r") as prediction:
        group = prediction["candidates"]
        theta_true = slopes_to_angles(source["truth/slope_xy_bins_per_z"][:])[0]
        xs = prediction["grid_x_start_bin"][:]
        ys = prediction["grid_y_start_bin"][:]
        size = int(prediction.attrs["crop_size"])
        selected = [
            i for i in range(len(group["event_id"]))
            if bool(group["truth_matched"][i])
            and theta_true[int(group["event_index"][i])] > 20.0
        ]
        selected.sort(key=lambda i: float(group["presence_score"][i]), reverse=True)
        for i in selected[:len(collect_background.cache)]:
            event = int(group["event_index"][i])
            x0 = int(xs[int(group["grid_col"][i])])
            y0 = int(ys[int(group["grid_row"][i])])
            raw = source["volumes_raw"][event, :, y0:y0 + size, x0:x0 + size]
            crops.append(normalize_per_slice(raw, source["background_mu"][event]))
    return np.stack(crops)


def score(model: torch.nn.Module, values: np.ndarray) -> np.ndarray:
    with torch.no_grad():
        logits = model(torch.from_numpy(np.ascontiguousarray(values[:, None])).float())["presence_logit"]
    return torch.sigmoid(logits).cpu().numpy()


def summarize(scores: dict[str, np.ndarray]) -> dict[str, object]:
    original = scores["original"]
    return {
        name: {
            "median": float(np.median(values)),
            "mean": float(np.mean(values)),
            "minimum": float(np.min(values)),
            "maximum": float(np.max(values)),
            "above_0p98": int((values >= 0.98).sum()),
            "median_change_from_original": float(np.median(values - original)),
        }
        for name, values in scores.items()
    }


def main() -> None:
    config = load_config(CONFIG)
    model = build_model(config).cpu()
    checkpoint = torch.load(CHECKPOINT, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    background = collect_background()
    collect_background.cache = background  # type: ignore[attr-defined]
    signal = collect_signal()
    rng = np.random.default_rng(20260827)
    permutation = rng.permutation(background.shape[1])

    result: dict[str, object] = {}
    for name, crops in (("background_candidates", background), ("matched_signal_theta_gt20", signal)):
        pixel_shuffled = crops.copy()
        for event in range(len(pixel_shuffled)):
            for z in range(pixel_shuffled.shape[1]):
                pixel_shuffled[event, z] = rng.permutation(pixel_shuffled[event, z].ravel()).reshape(20, 20)
        variants = {
            "original": crops,
            "z_reversed": crops[:, ::-1],
            "z_permuted": crops[:, permutation],
            "xy_pixels_permuted_per_layer": pixel_shuffled,
            "per_layer_spatial_mean_removed": (
                crops - crops.mean(axis=(2, 3), keepdims=True)
            ),
        }
        scores = {variant: score(model, values) for variant, values in variants.items()}
        result[name] = {"events": len(crops), "scores": summarize(scores)}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
