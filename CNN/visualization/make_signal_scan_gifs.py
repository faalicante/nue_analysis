#!/usr/bin/env python3
"""Create smoothed GIFs for missed and clearly detected signal scan events."""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from training.cnn21d.metrics import slopes_to_angles
from visualization.make_background_candidate_gifs import root_th2_smooth
from scanning.scan_cnn21d_volumes import cluster_contains_truth, connected_components


N_LAYERS = 57
PLATE_INDEX_OFFSET = 3
MAX_PROPAGATED_PLATES = 40
Z_BIN_SIZE_UM = 1350.0


def root_crop_start(edges: np.ndarray, center_um: float, size: int) -> int:
    """Match ROOT FindFixBin plus crop_bw's even-sized crop convention."""
    if size <= 0 or size % 2:
        raise ValueError("crop size must be a positive even integer")
    if center_um < edges[0] or center_um >= edges[-1]:
        raise ValueError(f"centre {center_um} outside histogram support")
    center_bin_one_based = int(np.searchsorted(edges, center_um, side="right"))
    first_bin_one_based = center_bin_one_based - size // 2
    start = first_bin_one_based - 1
    if start < 0 or start + size > len(edges) - 1:
        raise ValueError(f"crop {size} around {center_um} exceeds stored volume")
    return start


def mc_projection(source: h5py.File, event_index: int) -> dict[str, float | int]:
    """Recompute the crop centre using min(40, 57 - (truth_plate - 3))."""
    plate = int(source["truth/plate"][event_index])
    crop_plate = plate - PLATE_INDEX_OFFSET
    propagated_plates = min(MAX_PROPAGATED_PLATES, N_LAYERS - crop_plate)
    half_length_um = 0.5 * propagated_plates * Z_BIN_SIZE_UM
    center_x_um = (
        float(source["truth/xpos_um"][event_index])
        + float(source["truth/ltx_mrad_corrected"][event_index]) * half_length_um / 1000.0
    )
    center_y_um = (
        float(source["truth/ypos_um"][event_index])
        + float(source["truth/lty_mrad_corrected"][event_index]) * half_length_um / 1000.0
    )
    stored_x = float(source["truth/propagated_center_x_um"][event_index])
    stored_y = float(source["truth/propagated_center_y_um"][event_index])
    # The corrected slopes are stored as float32, whereas production computed
    # the centres from the original text precision before writing float64.
    if not np.allclose([center_x_um, center_y_um], [stored_x, stored_y], rtol=0, atol=1e-3):
        raise ValueError("recomputed MC projection differs from the stored production value")
    return {
        "truth_plate": plate,
        "crop_plate": crop_plate,
        "propagated_plates": propagated_plates,
        "center_x_um": stored_x,
        "center_y_um": stored_y,
        "recomputed_center_x_um_from_stored_float32_slope": center_x_um,
        "recomputed_center_y_um_from_stored_float32_slope": center_y_um,
    }


def render_frames(
    crop: np.ndarray,
    *,
    extent_mm: tuple[float, float, float, float],
    title: str,
) -> list[Image.Image]:
    finite = crop[np.isfinite(crop)]
    vmin = float(np.quantile(finite, 0.02))
    vmax = max(float(np.quantile(finite, 0.998)), vmin + 1.0)
    figure, axis = plt.subplots(figsize=(6.4, 6.1), constrained_layout=True)
    image = axis.imshow(
        crop[0], origin="lower", cmap="viridis", vmin=vmin, vmax=vmax,
        extent=extent_mm, interpolation="nearest", aspect="equal",
    )
    colorbar = figure.colorbar(image, ax=axis, shrink=0.84)
    colorbar.set_label("conteggi smussati")
    axis.set_xlabel("x [mm]")
    axis.set_ylabel("y [mm]")
    layer_title = axis.set_title(f"{title}\nlayer 01/57", fontsize=10.5)
    frames: list[Image.Image] = []
    for z in range(crop.shape[0]):
        image.set_data(crop[z])
        layer_title.set_text(f"{title}\nlayer {z + 1:02d}/57")
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=100)
        buffer.seek(0)
        frames.append(Image.open(buffer).convert("P", palette=Image.Palette.ADAPTIVE, colors=256))
        buffer.close()
    plt.close(figure)
    return frames


def save_gif(
    source: h5py.File,
    event_index: int,
    *,
    start_x: int,
    start_y: int,
    crop_size: int,
    title: str,
    output: Path,
    frame_duration_ms: int,
) -> dict[str, object]:
    raw = np.asarray(source["volumes_raw"][event_index])
    # Match the background GIFs: smooth the full 200x200 plane before cropping.
    smoothed = root_th2_smooth(raw)
    crop = smoothed[:, start_y:start_y + crop_size, start_x:start_x + crop_size]
    if crop.shape != (N_LAYERS, crop_size, crop_size):
        raise ValueError(f"unexpected crop shape {crop.shape}")
    x_edges = np.asarray(source["x_edges_um"][event_index])
    y_edges = np.asarray(source["y_edges_um"][event_index])
    extent_mm = (
        float(x_edges[start_x] / 1000.0),
        float(x_edges[start_x + crop_size] / 1000.0),
        float(y_edges[start_y] / 1000.0),
        float(y_edges[start_y + crop_size] / 1000.0),
    )
    frames = render_frames(crop, extent_mm=extent_mm, title=title)
    durations = [1000] + [frame_duration_ms] * (len(frames) - 2) + [1000]
    frames[0].save(
        output, save_all=True, append_images=frames[1:], duration=durations,
        loop=0, optimize=False, disposal=2,
    )
    return {
        "gif": str(output),
        "crop_size": crop_size,
        "start_x_bin_zero_based": start_x,
        "start_y_bin_zero_based": start_y,
        "extent_mm": list(extent_mm),
        "frames": len(frames),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--volumes", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.98)
    parser.add_argument("--theta-min", type=float, default=10.0)
    parser.add_argument("--clear-theta-min", type=float, default=20.0)
    parser.add_argument("--frame-duration-ms", type=int, default=160)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, object]] = []
    with h5py.File(args.volumes, "r") as source, h5py.File(args.predictions, "r") as prediction:
        scores = np.asarray(prediction["window_presence_score"])
        true_slopes = np.asarray(source["truth/slope_xy_bins_per_z"])
        theta_true = slopes_to_angles(true_slopes)[0]
        truth_x_bin = np.asarray(source["truth/propagated_center_x_bin"])
        truth_y_bin = np.asarray(source["truth/propagated_center_y_bin"])
        event_ids = np.asarray(source["event_id"])
        x_positions = np.asarray(prediction["grid_x_start_bin"], dtype=int)
        y_positions = np.asarray(prediction["grid_y_start_bin"], dtype=int)
        scan_crop_size = int(prediction.attrs["crop_size"])

        missed_indices: list[int] = []
        detected_indices: list[int] = []
        for event_index in range(len(event_ids)):
            components = connected_components(scores[event_index] >= args.threshold)
            detected = any(
                cluster_contains_truth(
                    component,
                    float(truth_x_bin[event_index]), float(truth_y_bin[event_index]),
                    x_positions.tolist(), y_positions.tolist(), scan_crop_size,
                )
                for component in components
            )
            if theta_true[event_index] > args.theta_min and not detected:
                missed_indices.append(event_index)
            if theta_true[event_index] > args.clear_theta_min and detected:
                detected_indices.append(event_index)

        if not missed_indices:
            raise ValueError(
                f"no missed theta>{args.theta_min:g} events found at threshold "
                f"{args.threshold:g}"
            )
        if not detected_indices:
            raise ValueError("no clearly detected signal event found")
        clear_index = max(detected_indices, key=lambda index: float(scores[index].max()))

        for kind, event_indices in (("clear", [clear_index]), ("missed", missed_indices)):
            event_indices = sorted(event_indices, key=lambda index: float(theta_true[index]), reverse=True)
            for rank, event_index in enumerate(event_indices, start=1):
                event_id = int(event_ids[event_index])
                projection = mc_projection(source, event_index)
                x_edges = np.asarray(source["x_edges_um"][event_index])
                y_edges = np.asarray(source["y_edges_um"][event_index])
                start_x = root_crop_start(x_edges, float(projection["center_x_um"]), 40)
                start_y = root_crop_start(y_edges, float(projection["center_y_um"]), 40)
                max_score = float(scores[event_index].max())
                label = "Signal evidente" if kind == "clear" else "Signal non rilevato"
                title = (
                    f"{label} | event_id={event_id} | theta_MC={theta_true[event_index]:.2f} mrad\n"
                    f"crop 40x40 su proiezione MC | min(40, 57-p0)="
                    f"{int(projection['propagated_plates'])} | max score={max_score:.6f}"
                )
                filename = (
                    f"{kind}_rank_{rank:02d}_event_{event_id:05d}_"
                    f"theta{theta_true[event_index]:06.2f}_mc40.gif"
                )
                output = args.output_dir / filename
                result = save_gif(
                    source, event_index, start_x=start_x, start_y=start_y,
                    crop_size=40, title=title, output=output,
                    frame_duration_ms=args.frame_duration_ms,
                )
                manifest_rows.append({
                    "kind": kind,
                    "rank": rank,
                    "event_index": int(event_index),
                    "event_id": event_id,
                    "theta_true_mrad": float(theta_true[event_index]),
                    "maximum_scan_score": max_score,
                    "mc_projection": projection,
                } | result)

                if kind == "clear":
                    grid_row, grid_col = np.unravel_index(
                        int(np.argmax(scores[event_index])), scores[event_index].shape
                    )
                    start_x_20 = int(x_positions[grid_col])
                    start_y_20 = int(y_positions[grid_row])
                    title20 = (
                        f"Signal evidente | event_id={event_id} | theta_MC={theta_true[event_index]:.2f} mrad\n"
                        f"crop CNN 20x20 a score massimo | P(signal)={max_score:.6f}"
                    )
                    output20 = args.output_dir / (
                        f"clear_event_{event_id:05d}_theta{theta_true[event_index]:06.2f}_maxscore20.gif"
                    )
                    result20 = save_gif(
                        source, event_index, start_x=start_x_20, start_y=start_y_20,
                        crop_size=20, title=title20, output=output20,
                        frame_duration_ms=args.frame_duration_ms,
                    )
                    manifest_rows.append({
                        "kind": "clear_max_score_20",
                        "event_index": int(event_index),
                        "event_id": event_id,
                        "theta_true_mrad": float(theta_true[event_index]),
                        "maximum_scan_score": max_score,
                        "grid_row": int(grid_row),
                        "grid_col": int(grid_col),
                        "mc_projection": projection,
                    } | result20)

    manifest = {
        "threshold": args.threshold,
        "missed_selection": f"theta_true > {args.theta_min:g} mrad and no truth-matched cluster",
        "clear_selection": (
            f"highest maximum scan score among truth-detected events with theta_true > "
            f"{args.clear_theta_min:g} mrad"
        ),
        "mc_center_formula": (
            "p0 = truth_plate - 3; n = min(40, 57 - p0); "
            "center_xy = vertex_xy + corrected_slope_mrad * 1350 um * 0.5 * n / 1000"
        ),
        "smoothing": "one ROOT TH2::Smooth(k5a)-equivalent pass on each complete 200x200 layer",
        "tagging_overlays": False,
        "frame_duration_ms": args.frame_duration_ms,
        "gifs": manifest_rows,
    }
    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "gif_count": len(manifest_rows),
        "missed_count": sum(row["kind"] == "missed" for row in manifest_rows),
        "output_dir": str(args.output_dir.resolve()),
        "manifest": str(manifest_path.resolve()),
    }, indent=2))


if __name__ == "__main__":
    main()
