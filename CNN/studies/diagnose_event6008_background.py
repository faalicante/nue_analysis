#!/usr/bin/env python3
"""Show why the old free fit for event 6008 followed background activity."""

from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import uproot


HDF5_PATH = Path("cnn_dataset.h5")
HDF5_INDEX = 9876
ROOT_PATH = Path("/Users/fabioali/cernbox/shift/nue_Euniform_FLUKA/b000021.0.0.6009.trk.root")
FIT_RESULTS = Path("runs/cnn21d_sample_type_no_poisson_theta10_50/pilot/centroid_fit/fit_results.json")
OUTPUT_DIR = Path("runs/cnn21d_sample_type_no_poisson_theta10_50/pilot/truth_tracks")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fit_records = json.loads(FIT_RESULTS.read_text())["records"]
    record = next(item for item in fit_records if item["hdf5_index"] == HDF5_INDEX)
    fit = record["raw32_fit"]

    with h5py.File(HDF5_PATH, "r") as hdf5:
        raw = np.asarray(hdf5["regions_raw"][HDF5_INDEX])
        mu = np.asarray(hdf5["background_mu"][HDF5_INDEX], dtype=np.float64)
        label = np.asarray(hdf5["slope_xy"][HDF5_INDEX], dtype=np.float64)
    significance = (raw - mu[:, None, None]) / np.sqrt(mu[:, None, None] + 1e-6)

    root_file = uproot.open(ROOT_PATH)
    histogram = root_file["XYPseg"]
    x_edges, y_edges, _ = (axis.edges() for axis in histogram.axes)

    # Same centre used for event 6008 in crop_bw.C.  Its remaining 34 plates
    # are below the new cap of 40, so the crop position is unchanged.
    crop_x_um = 364282.5902
    crop_y_um = 49679.4596
    center_x_bin = np.searchsorted(x_edges, crop_x_um, side="right")
    center_y_bin = np.searchsorted(y_edges, crop_y_um, side="right")
    first_x_bin = center_x_bin - 16
    first_y_bin = center_y_bin - 16

    names = ["s/s.eFlag", "s/s.eMCTrack", "s/s.eX", "s/s.eY", "s/s.eZ"]
    arrays = root_file["tracks"].arrays(names, library="ak")
    flag, mc_track, x_um, y_um, z_um = (
        ak.to_numpy(ak.flatten(arrays[name])) for name in names
    )
    x = np.searchsorted(x_edges, x_um, side="right") - first_x_bin
    y = np.searchsorted(y_edges, y_um, side="right") - first_y_bin
    z = np.rint((z_um + 75600.0) / 1350.0).astype(np.int64)
    inside = (x >= 0) & (x < 32) & (y >= 0) & (y < 32) & (z >= 0) & (z < 57)
    signal = inside & (flag == 1)
    background = inside & (flag == 0)
    primary = signal & (mc_track == 1)

    background_counts = np.zeros((57, 32, 32), dtype=np.int32)
    signal_counts = np.zeros_like(background_counts)
    np.add.at(background_counts, (z[background], y[background], x[background]), 1)
    np.add.at(signal_counts, (z[signal], y[signal], x[signal]), 1)

    zz = np.arange(57, dtype=np.float64)
    fit_x = fit["intercept_x"] + fit["slope_x"] * zz
    fit_y = fit["intercept_y"] + fit["slope_y"] * zz

    # Actual nominal electron line, anchored at the MC interaction vertex.
    interaction_layer = 23
    vertex_x_bin = np.searchsorted(x_edges, 364657.26, side="right") - first_x_bin
    vertex_y_bin = np.searchsorted(y_edges, 50673.81, side="right") - first_y_bin
    truth_x = vertex_x_bin + label[0] * (zz - interaction_layer)
    truth_y = vertex_y_bin + label[1] * (zz - interaction_layer)

    figure, axes = plt.subplots(3, 2, figsize=(14, 13), constrained_layout=True, sharex=True, sharey="row")
    total_projections = (significance.max(axis=1).T, significance.max(axis=2).T)
    background_projections = (background_counts.max(axis=1).T, background_counts.max(axis=2).T)
    signal_projections = (signal_counts.max(axis=1).T, signal_counts.max(axis=2).T)
    coordinates = ((x, fit_x, truth_x), (y, fit_y, truth_y))
    labels = ("x", "y")

    for column, (coord, old_line, truth_line) in enumerate(coordinates):
        axis = axes[0, column]
        axis.imshow(total_projections[column], origin="lower", aspect="auto", cmap="magma", extent=(0, 57, 0, 32))
        axis.plot(zz, old_line, color="lime", lw=2.5, label="vecchio fit RANSAC")
        axis.plot(zz, truth_line, "--", color="cyan", lw=2.5, label="asse vero ancorato al vertice")
        axis.scatter(z[signal], coord[signal], s=14, facecolors="none", edgecolors="white", alpha=0.55, linewidths=0.7, label="tutto eFlag=1")
        axis.scatter(z[primary], coord[primary], s=55, c="red", edgecolors="white", linewidths=0.7, label="elettrone eMCTrack=1")
        axis.set_title(f"HDF5 totale — proiezione z–{labels[column]}")
        axis.legend(fontsize=8, loc="best")

        axis = axes[1, column]
        projection = background_projections[column]
        axis.imshow(projection, origin="lower", aspect="auto", cmap="magma", norm=LogNorm(vmin=1, vmax=max(2, int(projection.max()))), extent=(0, 57, 0, 32))
        axis.plot(zz, old_line, color="lime", lw=2.5)
        axis.set_title(f"Solo fondo: s.eFlag=0 — z–{labels[column]}")

        axis = axes[2, column]
        projection = signal_projections[column]
        axis.imshow(projection, origin="lower", aspect="auto", cmap="Blues", vmin=0, vmax=max(1, int(projection.max())), extent=(0, 57, 0, 32))
        axis.scatter(z[signal], coord[signal], s=12, c="royalblue", alpha=0.45, linewidths=0)
        axis.scatter(z[primary], coord[primary], s=55, c="red", edgecolors="black", linewidths=0.5)
        axis.plot(zz, truth_line, "--", color="cyan", lw=2.5)
        axis.set_title(f"Solo segnale: s.eFlag=1 — z–{labels[column]}")
        axis.set_xlabel("layer z")

    for row in range(3):
        axes[row, 0].set_ylabel("coordinata x [bin]")
        axes[row, 1].set_ylabel("coordinata y [bin]")
        for column in range(2):
            axes[row, column].set_xlim(0, 56)
            axes[row, column].set_ylim(0, 32)

    # Quantify whether the significant voxels selected by the old line contain
    # any true signal segment.
    threshold_z, threshold_y, threshold_x = np.where(significance > 5.0)
    distance = np.hypot(
        threshold_x - (fit["intercept_x"] + fit["slope_x"] * threshold_z),
        threshold_y - (fit["intercept_y"] + fit["slope_y"] * threshold_z),
    )
    old_inliers = distance < 1.6
    signal_voxels = set(zip(z[signal].tolist(), y[signal].tolist(), x[signal].tolist()))
    primary_voxels = set(zip(z[primary].tolist(), y[primary].tolist(), x[primary].tolist()))
    old_voxels = list(zip(threshold_z[old_inliers], threshold_y[old_inliers], threshold_x[old_inliers]))
    old_signal_overlap = sum(tuple(map(int, voxel)) in signal_voxels for voxel in old_voxels)
    old_primary_overlap = sum(tuple(map(int, voxel)) in primary_voxels for voxel in old_voxels)

    figure.suptitle(
        "Evento 6008 / HDF5 9876 — il vecchio fit segue una struttura di fondo\n"
        f"inlier >5σ del fit: {len(old_voxels)}; sovrapposti a eFlag=1: {old_signal_overlap}; "
        f"sovrapposti al primario: {old_primary_overlap}",
        fontsize=15,
    )
    output_plot = OUTPUT_DIR / "event_6008_background_vs_electron_overlay.png"
    figure.savefig(output_plot, dpi=180)
    plt.close(figure)

    result = {
        "hdf5_index": HDF5_INDEX,
        "root_file": str(ROOT_PATH),
        "segments_inside_crop": {
            "all": int(inside.sum()),
            "background_eflag0": int(background.sum()),
            "signal_eflag1": int(signal.sum()),
            "primary_electron": int(primary.sum()),
        },
        "old_fit": {
            "slope_bins_per_layer": [fit["slope_x"], fit["slope_y"]],
            "theta_mrad": fit["theta_mrad"],
            "significant_inlier_voxels": len(old_voxels),
            "inlier_voxels_overlapping_signal": old_signal_overlap,
            "inlier_voxels_overlapping_primary": old_primary_overlap,
        },
        "label_slope_bins_per_layer": label.tolist(),
        "primary_points_inside_crop": list(zip(z[primary].tolist(), x[primary].tolist(), y[primary].tolist())),
        "plot": str(output_plot),
    }
    output_json = OUTPUT_DIR / "event_6008_background_vs_electron.json"
    output_json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
