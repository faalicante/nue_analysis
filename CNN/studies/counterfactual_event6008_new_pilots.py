#!/usr/bin/env python3
"""Remove eFlag==1 tracks from event 6008 and re-evaluate both new pilots."""

from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import h5py
import numpy as np
import torch
import uproot

from training.cnn_dataset import extract_xy_crop, normalize_per_slice
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, load_config


HDF5_PATH = Path("cnn_dataset_signal_t5_center40.h5")
EVENT_ID = 6008
ROOT_PATH = Path("/Users/fabioali/cernbox/shift/nue_Euniform_FLUKA/b000021.0.0.6009.trk.root")
OUTPUT_PATH = Path("runs/cnn21d_signal_t5_center40/event6008_counterfactual.json")


def infer(model: torch.nn.Module, raw: np.ndarray, mu: float, crop_size: int) -> dict[str, object]:
    crop = extract_xy_crop(raw, crop_size=crop_size)
    volume = torch.from_numpy(normalize_per_slice(crop, mu)[None, None])
    with torch.no_grad():
        output = model(volume)
    probability = float(torch.sigmoid(output["presence_logit"])[0])
    slope = output["slope_xy"][0].numpy()
    theta = float(slopes_to_angles(slope[None])[0][0])
    return {"probability_signal": probability, "predicted_signal": probability >= 0.5,
            "slope_pred": slope.tolist(), "theta_pred_mrad": theta}


def main() -> None:
    with h5py.File(HDF5_PATH, "r") as hdf5:
        ids = np.asarray(hdf5["signal_event_id"][:], dtype=np.int64)
        types = np.asarray(hdf5["sample_type"][:], dtype=np.int64)
        index = int(np.flatnonzero((ids == EVENT_ID) & (types == 2))[0])
        raw = np.asarray(hdf5["regions_raw"][index], dtype=np.float32)
        mu = float(np.asarray(hdf5["background_mu"][index]).reshape(-1)[0])

    root_file = uproot.open(ROOT_PATH)
    histogram = root_file["XYPseg"]
    x_edges, y_edges, _ = (axis.edges() for axis in histogram.axes)
    crop_x_um, crop_y_um = 364282.5902, 49679.4596
    first_x_bin = np.searchsorted(x_edges, crop_x_um, side="right") - 16
    first_y_bin = np.searchsorted(y_edges, crop_y_um, side="right") - 16
    names = ["s/s.eFlag", "s/s.eX", "s/s.eY", "s/s.eZ"]
    arrays = root_file["tracks"].arrays(names, library="ak")
    flag, x_um, y_um, z_um = (ak.to_numpy(ak.flatten(arrays[name])) for name in names)
    x = np.searchsorted(x_edges, x_um, side="right") - first_x_bin
    y = np.searchsorted(y_edges, y_um, side="right") - first_y_bin
    z = np.rint((z_um + 75600.0) / 1350.0).astype(np.int64)
    signal = (flag == 1) & (x >= 0) & (x < 32) & (y >= 0) & (y < 32) & (z >= 0) & (z < 57)
    signal_counts = np.zeros_like(raw, dtype=np.float32)
    np.add.at(signal_counts, (z[signal], y[signal], x[signal]), 1)
    residual = raw - signal_counts
    negative_voxels = int((residual < 0).sum())
    background_only = np.maximum(residual, 0)

    result: dict[str, object] = {
        "event_id": EVENT_ID, "hdf5_index": index,
        "eflag1_segments_inside_crop": int(signal.sum()),
        "signal_counts_sum": float(signal_counts.sum()),
        "negative_voxels_before_clip": negative_voxels,
        "pilots": {},
    }
    for name, config_path, checkpoint_path in (
        ("crop20_all_t5", Path("training/configs/cnn21d_signal_t5_center40_all_crop20.yaml"), Path("runs/cnn21d_signal_t5_center40/all_t5/crop20/pilot/best_model.pt")),
        ("crop32_all_t5", Path("training/configs/cnn21d_signal_t5_center40_all_crop32.yaml"), Path("runs/cnn21d_signal_t5_center40/all_t5/crop32/pilot/best_model.pt")),
    ):
        config = load_config(config_path)
        model = build_model(config).cpu()
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
        model.load_state_dict(checkpoint["model_state"]); model.eval()
        crop_size = int(config["data"]["dataset"]["crop_size"])
        full = infer(model, raw, mu, crop_size)
        removed = infer(model, background_only, mu, crop_size)
        result["pilots"][name] = {
            "full_event": full,
            "eflag1_removed": removed,
            "delta_probability_full_minus_removed": float(full["probability_signal"]) - float(removed["probability_signal"]),
        }
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
