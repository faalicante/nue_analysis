#!/usr/bin/env python3
"""Inspect signal-only FEDRA segments for selected nue events.

The ROOT selection is strictly ``s.eFlag == 1``.  The primary electron
(``s.eMCTrack == 1``) is highlighted separately from all other particles of
the same MC interaction.  Plots are expressed in the same 50 um xy bins and
57-layer convention used by the CNN crops.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.patches import Rectangle


DZ_UM = 1350.0
XY_BIN_UM = 50.0
N_LAYERS = 57
MAX_PROPAGATED_PLATES = 40
DZ_SHRINK = 1315.0 / 1350.0

BRANCHES = {
    "flag": "s/s.eFlag",
    "mc_event": "s/s.eMCEvt",
    "mc_track": "s/s.eMCTrack",
    "pid": "s/s.ePID",
    "x": "s/s.eX",
    "y": "s/s.eY",
    "z": "s/s.eZ",
    "tx": "s/s.eTX",
    "ty": "s/s.eTY",
}


def read_event_metadata(path: Path, requested: set[int]) -> dict[int, dict[str, float | int]]:
    records: dict[int, dict[str, float | int]] = {}
    with path.open(newline="") as stream:
        for fields in csv.reader(stream):
            event = int(fields[0])
            if event not in requested:
                continue
            records[event] = {
                "event": event,
                "vertex_x_um": float(fields[1]),
                "vertex_y_um": float(fields[2]),
                "tx_mrad_raw": float(fields[5]),
                "ty_mrad_raw": float(fields[6]),
                "energy_gev": float(fields[7]),
                "interaction_plate_file": int(fields[8]),
                "interaction_layer": int(fields[8]) - 3,
            }
    missing = requested.difference(records)
    if missing:
        raise RuntimeError(f"Events missing from {path}: {sorted(missing)}")
    return records


def flatten_selected(tree: uproot.TTree) -> tuple[dict[str, np.ndarray], dict[int, int]]:
    arrays = tree.arrays(list(BRANCHES.values()), library="ak")
    flags = ak.to_numpy(ak.flatten(arrays[BRANCHES["flag"]]))
    values, counts = np.unique(flags, return_counts=True)
    flag_counts = {int(value): int(count) for value, count in zip(values, counts, strict=True)}
    mask = arrays[BRANCHES["flag"]] == 1
    selected = {
        name: ak.to_numpy(ak.flatten(arrays[branch][mask]))
        for name, branch in BRANCHES.items()
        if name != "flag"
    }
    return selected, flag_counts


def angle_mrad(tx: float, ty: float) -> float:
    return float(1000.0 * np.arctan(np.hypot(tx, ty)))


def line_fit(z: np.ndarray, coordinate: np.ndarray) -> tuple[float, float]:
    slope, intercept = np.polyfit(z.astype(np.float64), coordinate.astype(np.float64), 1)
    return float(slope), float(intercept)


def analyze_event(
    root_path: Path,
    metadata: dict[str, float | int],
    output_dir: Path,
) -> dict[str, object]:
    event = int(metadata["event"])
    root_file = uproot.open(root_path)
    tree = root_file["tracks"]  # newest cycle
    data, flag_counts = flatten_selected(tree)

    mc_events, mc_event_counts = np.unique(data["mc_event"], return_counts=True)
    if mc_events.tolist() != [event]:
        raise RuntimeError(f"{root_path}: eFlag==1 has MC events {mc_events.tolist()}, expected {event}")

    primary_mask = data["mc_track"] == 1
    if primary_mask.sum() < 3:
        raise RuntimeError(f"{root_path}: fewer than three eFlag==1, eMCTrack==1 segments")
    fit_tx, fit_x0 = line_fit(data["z"][primary_mask], data["x"][primary_mask])
    fit_ty, fit_y0 = line_fit(data["z"][primary_mask], data["y"][primary_mask])

    tx_raw = float(metadata["tx_mrad_raw"]) / 1000.0
    ty_raw = float(metadata["ty_mrad_raw"]) / 1000.0
    tx_label = tx_raw * DZ_SHRINK
    ty_label = ty_raw * DZ_SHRINK
    layer0 = int(metadata["interaction_layer"])
    vertex_x = float(metadata["vertex_x_um"])
    vertex_y = float(metadata["vertex_y_um"])

    # This reproduces crop_bw.C exactly: propagate across at most 40 plates,
    # then place the fixed crop at the midpoint (at most 20 plates forward).
    propagated_plates = min(MAX_PROPAGATED_PLATES, N_LAYERS - layer0)
    crop_x = vertex_x + tx_label * DZ_UM * 0.5 * propagated_plates
    crop_y = vertex_y + ty_label * DZ_UM * 0.5 * propagated_plates

    # FEDRA stores the 57 plates at z=-75600,...,0 um.  This must not be
    # inferred from the first signal hit: late interactions (e.g. event 6008)
    # have no signal in the upstream plates.
    stack_z_min = -DZ_UM * (N_LAYERS - 1)
    layer = (data["z"] - stack_z_min) / DZ_UM
    x_bin = (data["x"] - crop_x) / XY_BIN_UM
    y_bin = (data["y"] - crop_y) / XY_BIN_UM
    nominal_x_bin = (
        vertex_x + tx_raw * DZ_UM * (layer - layer0) - crop_x
    ) / XY_BIN_UM
    nominal_y_bin = (
        vertex_y + ty_raw * DZ_UM * (layer - layer0) - crop_y
    ) / XY_BIN_UM

    in_crop32 = (np.abs(x_bin) < 16) & (np.abs(y_bin) < 16)
    in_crop20 = (np.abs(x_bin) < 10) & (np.abs(y_bin) < 10)
    layers = np.arange(N_LAYERS)
    count32 = np.asarray([np.sum(in_crop32 & np.isclose(layer, value)) for value in layers])
    count20 = np.asarray([np.sum(in_crop20 & np.isclose(layer, value)) for value in layers])

    figure, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    color = layer
    common = dict(c=color, cmap="viridis", s=11, alpha=0.48, linewidths=0)

    axes[0, 0].scatter(layer, x_bin, **common)
    axes[0, 0].scatter(layer[primary_mask], x_bin[primary_mask], s=34, c="crimson", label="eMCTrack=1")
    order = np.argsort(layer)
    axes[0, 0].plot(layer[order], nominal_x_bin[order], "--", color="deepskyblue", lw=2, label="direzione nominale")
    axes[0, 0].axhspan(-16, 16, color="gold", alpha=0.08, label="crop 32")
    axes[0, 0].axhspan(-10, 10, color="limegreen", alpha=0.10, label="crop CNN 20")
    axes[0, 0].set(xlabel="layer (0–56)", ylabel="x − centro crop [bin da 50 µm]", title="Proiezione x–z")
    axes[0, 0].legend(fontsize=8, loc="best")

    axes[0, 1].scatter(layer, y_bin, **common)
    axes[0, 1].scatter(layer[primary_mask], y_bin[primary_mask], s=34, c="crimson")
    axes[0, 1].plot(layer[order], nominal_y_bin[order], "--", color="deepskyblue", lw=2)
    axes[0, 1].axhspan(-16, 16, color="gold", alpha=0.08)
    axes[0, 1].axhspan(-10, 10, color="limegreen", alpha=0.10)
    axes[0, 1].set(xlabel="layer (0–56)", ylabel="y − centro crop [bin da 50 µm]", title="Proiezione y–z")

    points = axes[1, 0].scatter(x_bin, y_bin, **common)
    axes[1, 0].scatter(x_bin[primary_mask], y_bin[primary_mask], s=38, c="crimson", label="elettrone primario")
    axes[1, 0].add_patch(Rectangle((-16, -16), 32, 32, fill=False, ec="goldenrod", lw=2, label="32×32"))
    axes[1, 0].add_patch(Rectangle((-10, -10), 20, 20, fill=False, ec="limegreen", lw=2, label="20×20"))
    axes[1, 0].set_aspect("equal", adjustable="box")
    axes[1, 0].set(xlabel="x − centro crop [bin]", ylabel="y − centro crop [bin]", title="Vista xy: tutti i segmenti signal-only")
    axes[1, 0].legend(fontsize=8, loc="best")
    figure.colorbar(points, ax=axes[1, 0], label="layer")

    axes[1, 1].step(layers, count32, where="mid", color="goldenrod", lw=2, label="dentro 32×32")
    axes[1, 1].step(layers, count20, where="mid", color="limegreen", lw=2, label="dentro 20×20")
    axes[1, 1].set(xlabel="layer (0–56)", ylabel="segmenti con eFlag=1", title="Segnale realmente visibile nel crop")
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.25)

    figure.suptitle(
        f"Evento {event} — solo s.eFlag==1 | primario: "
        f"({1000*fit_tx:+.1f}, {1000*fit_ty:+.1f}) mrad, θ={angle_mrad(fit_tx, fit_ty):.1f} mrad",
        fontsize=14,
    )
    plot_path = output_dir / f"event_{event}_signal_truth_tracks.png"
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)

    primary_tx_median = float(np.median(data["tx"][primary_mask]))
    primary_ty_median = float(np.median(data["ty"][primary_mask]))
    result: dict[str, object] = {
        "event": event,
        "root_file": str(root_path),
        "tree_cycle": tree.object_path,
        "tree_entries": int(tree.num_entries),
        "all_segment_flag_counts": flag_counts,
        "signal_selection": "s.eFlag == 1",
        "signal_segments": int(data["x"].size),
        "signal_mc_event_counts": {str(int(key)): int(value) for key, value in zip(mc_events, mc_event_counts, strict=True)},
        "signal_mc_tracks": int(np.unique(data["mc_track"]).size),
        "signal_layers": int(np.unique(data["pid"]).size),
        "primary_electron": {
            "selection": "s.eFlag == 1 && s.eMCTrack == 1",
            "segments": int(primary_mask.sum()),
            "fit_tx_mrad": 1000.0 * fit_tx,
            "fit_ty_mrad": 1000.0 * fit_ty,
            "median_segment_tx_mrad": 1000.0 * primary_tx_median,
            "median_segment_ty_mrad": 1000.0 * primary_ty_median,
            "fit_theta_mrad": angle_mrad(fit_tx, fit_ty),
            "fit_slope_bins_per_layer": [DZ_UM * fit_tx / XY_BIN_UM, DZ_UM * fit_ty / XY_BIN_UM],
        },
        "metadata_direction": {
            "raw_tx_mrad": float(metadata["tx_mrad_raw"]),
            "raw_ty_mrad": float(metadata["ty_mrad_raw"]),
            "raw_theta_mrad": angle_mrad(tx_raw, ty_raw),
            "label_tx_mrad_after_1315_over_1350": 1000.0 * tx_label,
            "label_ty_mrad_after_1315_over_1350": 1000.0 * ty_label,
            "label_theta_mrad": angle_mrad(tx_label, ty_label),
            "hdf5_slope_bins_per_layer": [DZ_UM * tx_label / XY_BIN_UM, DZ_UM * ty_label / XY_BIN_UM],
        },
        "crop": {
            "max_propagated_plates": MAX_PROPAGATED_PLATES,
            "propagated_plates": propagated_plates,
            "center_offset_plates": 0.5 * propagated_plates,
            "center_x_um": crop_x,
            "center_y_um": crop_y,
            "signal_segments_in_32": int(in_crop32.sum()),
            "signal_layers_in_32": int(np.count_nonzero(count32)),
            "signal_segments_in_20": int(in_crop20.sum()),
            "signal_layers_in_20": int(np.count_nonzero(count20)),
            "nominal_full_stack_displacement_bins": [
                DZ_UM * tx_raw * (N_LAYERS - 1) / XY_BIN_UM,
                DZ_UM * ty_raw * (N_LAYERS - 1) / XY_BIN_UM,
            ],
        },
        "plot": str(plot_path),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("/Users/fabioali/cernbox/shift/nue_Euniform_FLUKA"))
    parser.add_argument("--events", type=int, nargs="+", default=[938, 6008])
    parser.add_argument("--output-dir", type=Path, default=Path("runs/cnn21d_sample_type_no_poisson_theta10_50/pilot/truth_tracks"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    requested = set(args.events)
    metadata = read_event_metadata(args.input_dir / "nue_int_10k.txt", requested)
    results = []
    for event in args.events:
        root_path = args.input_dir / f"b000021.0.0.{event + 1}.trk.root"
        results.append(analyze_event(root_path, metadata[event], args.output_dir))

    output = args.output_dir / "truth_track_results.json"
    output.write_text(json.dumps({"events": results}, indent=2) + "\n")
    print(json.dumps({"output": str(output), "events": results}, indent=2))


if __name__ == "__main__":
    main()
