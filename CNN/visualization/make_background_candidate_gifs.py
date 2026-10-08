#!/usr/bin/env python3
"""Create one smoothed 57-layer GIF for each background scan candidate."""

from __future__ import annotations

import argparse
import csv
import glob
import io
import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np
from PIL import Image
from scipy.ndimage import convolve, label
from scipy.stats import poisson

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT_K5A = np.asarray(
    [
        [0, 0, 1, 0, 0],
        [0, 2, 2, 2, 0],
        [1, 2, 5, 2, 1],
        [0, 2, 2, 2, 0],
        [0, 0, 1, 0, 0],
    ],
    dtype=np.float32,
)

POISSON_TAIL_PROBABILITY = 1.0e-4
CNN_RIDGE_RADIUS_BINS = 1.0
CNN_RIDGE_MIN_LAYERS = 4


def root_th2_smooth(layers: np.ndarray) -> np.ndarray:
    """Vectorized equivalent of one ROOT TH2::Smooth("k5a") pass."""
    values = np.asarray(layers, dtype=np.float32)
    kernel = ROOT_K5A[None]
    numerator = convolve(values, kernel, mode="constant", cval=0.0)
    denominator = convolve(np.ones_like(values), kernel, mode="constant", cval=0.0)
    return numerator / denominator


def centered_window_start(
    inner_start: int,
    inner_size: int,
    display_size: int,
    full_size: int,
) -> int:
    """Center a display window on a CNN crop, clamping it at cell boundaries."""
    if not 0 < inner_size <= display_size <= full_size:
        raise ValueError("expected 0 < CNN crop <= display crop <= full cell")
    center = inner_start + inner_size / 2.0
    requested = int(round(center - display_size / 2.0))
    return min(max(requested, 0), full_size - display_size)


def fitted_display_crop_size(
    requested_size: int,
    cnn_crop_size: int,
    height: int,
    width: int,
) -> int:
    """Fit a square display crop inside a possibly clipped scan cell."""
    if requested_size < cnn_crop_size:
        raise ValueError("display crop must be at least as large as the CNN crop")
    available = min(height, width)
    if available < cnn_crop_size:
        raise ValueError("scan cell is smaller than the CNN crop")
    return min(requested_size, available)


def component_at_score_threshold(
    scores: np.ndarray, threshold: float, row: int, col: int,
) -> set[tuple[int, int]]:
    """Return the eight-connected score component containing one grid point."""
    labels, _ = label(np.asarray(scores) >= threshold, structure=np.ones((3, 3)))
    component_id = int(labels[row, col])
    if component_id == 0:
        return set()
    return {tuple(point) for point in np.argwhere(labels == component_id)}


def components_by_event(
    prediction_glob: str, threshold: float,
) -> dict[int, list[set[tuple[int, int]]]]:
    """Load score components used to exclude candidates found by another model."""
    result: dict[int, list[set[tuple[int, int]]]] = {}
    for path in sorted(map(Path, glob.glob(prediction_glob))):
        with h5py.File(path, "r") as prediction:
            for event_id, scores in zip(
                prediction["event_id"][:], prediction["window_presence_score"][:], strict=True,
            ):
                key = int(event_id)
                if key in result:
                    raise ValueError(f"duplicate event ID {key} in overlap predictions")
                labels, count = label(
                    np.asarray(scores) >= threshold, structure=np.ones((3, 3)),
                )
                result[key] = [
                    {tuple(point) for point in np.argwhere(labels == component_id)}
                    for component_id in range(1, count + 1)
                ]
    if not result:
        raise ValueError("overlap prediction glob is empty")
    return result


def poisson_hot_threshold(mu: float, alpha: float = POISSON_TAIL_PROBABILITY) -> int:
    """Smallest integer k such that P(Poisson(mu) >= k) <= alpha."""
    k = max(1, int(poisson.isf(alpha, mu)) + 1)
    while poisson.sf(k - 1, mu) > alpha:
        k += 1
    while k > 1 and poisson.sf(k - 2, mu) <= alpha:
        k -= 1
    return k


def cnn_consistent_hot_ridge(
    raw: np.ndarray,
    background_mu: float,
    predicted_slope_xy: np.ndarray,
) -> dict[str, object] | None:
    """Find a q-hot ridge following the CNN direction, even across z gaps."""
    values = np.asarray(raw, dtype=np.float64)
    slope = np.asarray(predicted_slope_xy, dtype=np.float64).reshape(2)
    threshold = poisson_hot_threshold(background_mu)
    z_index, y_index, x_index = np.nonzero(values >= threshold)
    if z_index.size < CNN_RIDGE_MIN_LAYERS:
        return None

    weights = np.maximum(values[z_index, y_index, x_index] - background_mu, 0.0)
    intercept_x = x_index - slope[0] * z_index
    intercept_y = y_index - slope[1] * z_index
    best_mask: np.ndarray | None = None
    best_key: tuple[int, float, int] | None = None
    for start_x, start_y in zip(intercept_x, intercept_y, strict=True):
        residual = np.hypot(
            x_index - (start_x + slope[0] * z_index),
            y_index - (start_y + slope[1] * z_index),
        )
        ridge_mask = residual <= CNN_RIDGE_RADIUS_BINS
        active_layers = np.unique(z_index[ridge_mask])
        key = (
            int(active_layers.size),
            float(weights[ridge_mask].sum()),
            int(ridge_mask.sum()),
        )
        if best_key is None or key > best_key:
            best_key = key
            best_mask = ridge_mask
    if best_mask is None or best_key is None or best_key[0] < CNN_RIDGE_MIN_LAYERS:
        return None

    mask = np.zeros(values.shape, dtype=bool)
    mask[z_index[best_mask], y_index[best_mask], x_index[best_mask]] = True
    active = np.flatnonzero(np.any(mask, axis=(1, 2)))
    residual = np.maximum(values, 0.0) * mask
    layer_y, layer_x = np.indices(values.shape[1:])
    centroid_x: list[float] = []
    centroid_y: list[float] = []
    for z in active:
        layer_weights = residual[z]
        layer_sum = float(layer_weights.sum())
        centroid_x.append(float((layer_weights * layer_x).sum() / layer_sum))
        centroid_y.append(float((layer_weights * layer_y).sum() / layer_sum))
    fitted_slope = np.asarray([
        np.polyfit(active, centroid_x, 1)[0],
        np.polyfit(active, centroid_y, 1)[0],
    ])
    return {
        "mask": mask,
        "charge": float((values - background_mu)[mask].sum()),
        "voxels": int(mask.sum()),
        "layers": int(active.size),
        "slope": fitted_slope,
    }


def estimate_shower_start(
    raw_crop: np.ndarray,
    background_mu: float,
    x_centers_um: np.ndarray,
    y_centers_um: np.ndarray,
    *,
    predicted_slope_xy: np.ndarray | None = None,
    candidate_center_xy_um: np.ndarray | None = None,
) -> dict[str, object]:
    """Estimate the first layer/centroid of the candidate-compatible 3D component."""
    values = np.asarray(raw_crop, dtype=np.float64)
    if values.ndim != 3:
        raise ValueError("raw crop must have [z,y,x] shape")
    if x_centers_um.shape != (values.shape[2],) or y_centers_um.shape != (values.shape[1],):
        raise ValueError("coordinate center arrays disagree with raw crop")
    threshold = poisson_hot_threshold(background_mu)
    hot = values >= threshold
    components, count = label(hot, structure=np.ones((3, 3, 3), dtype=bool))
    residual = values - background_mu
    candidates: list[dict[str, object]] = []
    y_indices, x_indices = np.indices(values.shape[1:])
    for component_id in range(1, count + 1):
        mask = components == component_id
        charge = float(residual[mask].sum())
        active = np.flatnonzero(np.any(mask, axis=(1, 2)))
        weights_3d = np.maximum(residual, 0.0) * mask
        weight_sum_3d = float(weights_3d.sum())
        component_x_um = float(
            (weights_3d * x_centers_um[None, None, :]).sum() / weight_sum_3d
        )
        component_y_um = float(
            (weights_3d * y_centers_um[None, :, None]).sum() / weight_sum_3d
        )
        fitted_slope = np.asarray([np.nan, np.nan], dtype=np.float64)
        if active.size >= 3:
            centroid_x: list[float] = []
            centroid_y: list[float] = []
            for z_index in active:
                layer_weights = np.maximum(residual[z_index], 0.0) * mask[z_index]
                layer_sum = float(layer_weights.sum())
                centroid_x.append(float((layer_weights * x_indices).sum() / layer_sum))
                centroid_y.append(float((layer_weights * y_indices).sum() / layer_sum))
            fitted_slope = np.asarray([
                np.polyfit(active, centroid_x, 1)[0],
                np.polyfit(active, centroid_y, 1)[0],
            ])
        candidates.append({
            "mask": mask,
            "charge": charge,
            "voxels": int(mask.sum()),
            "layers": int(active.size),
            "slope": fitted_slope,
            "component_x_um": component_x_um,
            "component_y_um": component_y_um,
        })
    if not candidates:
        return {
            "start_layer": None,
            "start_root_z_bin": None,
            "start_x_um": None,
            "start_y_um": None,
            "start_component_voxels": 0,
            "start_component_layers": 0,
            "start_component_charge": 0.0,
            "start_quality": "not_found",
            "start_at_first_observed_layer": False,
            "poisson_hot_count_threshold": threshold,
            "component_selection": "not_found",
        }
    selection = "maximum_charge"
    selected = max(candidates, key=lambda row: float(row["charge"]))
    if predicted_slope_xy is not None and candidate_center_xy_um is not None:
        ridge = cnn_consistent_hot_ridge(values, background_mu, predicted_slope_xy)
        if ridge is not None:
            selected = ridge
            selection = "cnn_consistent_qhot_ridge"
        else:
            predicted_slope = np.asarray(predicted_slope_xy, dtype=np.float64).reshape(2)
            candidate_center = np.asarray(candidate_center_xy_um, dtype=np.float64).reshape(2)
            bin_width_um = float(np.mean([
                np.median(np.diff(x_centers_um)), np.median(np.diff(y_centers_um))
            ]))
            eligible = [
                row for row in candidates
                if int(row["voxels"]) >= 8
                and int(row["layers"]) >= 3
                and np.all(np.isfinite(row["slope"]))
            ]
            if eligible:
                def association_cost(row: dict[str, object]) -> float:
                    slope_distance = float(np.linalg.norm(row["slope"] - predicted_slope))
                    center_distance_bins = float(np.hypot(
                        float(row["component_x_um"]) - candidate_center[0],
                        float(row["component_y_um"]) - candidate_center[1],
                    ) / bin_width_um)
                    return slope_distance + 0.01 * center_distance_bins

                selected = min(eligible, key=association_cost)
                selection = "predicted_slope_plus_candidate_center"
    best_mask = np.asarray(selected["mask"], dtype=bool)
    best_charge = float(selected["charge"])
    active_layers = np.flatnonzero(np.any(best_mask, axis=(1, 2)))
    start_index = int(active_layers[0])
    start_mask = best_mask[start_index]
    weights = np.maximum(residual[start_index], 0.0) * start_mask
    weight_sum = float(weights.sum())
    start_x_um = float((weights * x_centers_um[None, :]).sum() / weight_sum)
    start_y_um = float((weights * y_centers_um[:, None]).sum() / weight_sum)
    voxels = int(best_mask.sum())
    layers = int(active_layers.size)
    return {
        # GIF layer 1 corresponds to ROOT XYPseg z bin 2.
        "start_layer": start_index + 1,
        "start_root_z_bin": start_index + 2,
        "start_x_um": start_x_um,
        "start_y_um": start_y_um,
        "start_component_voxels": voxels,
        "start_component_layers": layers,
        "start_component_charge": best_charge,
        "start_quality": "coherent" if voxels >= 8 and layers >= 3 else "weak",
        "start_at_first_observed_layer": start_index == 0,
        "poisson_hot_count_threshold": threshold,
        "component_selection": selection,
        "selected_component_slope_x": float(selected["slope"][0]),
        "selected_component_slope_y": float(selected["slope"][1]),
    }


def render_frames(
    smoothed_crop: np.ndarray,
    *,
    extent_mm: tuple[float, float, float, float],
    title: str,
) -> list[Image.Image]:
    finite = smoothed_crop[np.isfinite(smoothed_crop)]
    vmin = float(np.quantile(finite, 0.02))
    vmax = max(float(np.quantile(finite, 0.998)), vmin + 1.0)
    figure, axis = plt.subplots(figsize=(6.4, 6.1), constrained_layout=True)
    image = axis.imshow(
        smoothed_crop[0], origin="lower", cmap="viridis",
        vmin=vmin, vmax=vmax, extent=extent_mm,
        interpolation="nearest", aspect="equal",
    )
    colorbar = figure.colorbar(image, ax=axis, shrink=0.84)
    colorbar.set_label("conteggi smussati")
    axis.set_xlabel("x [mm]")
    axis.set_ylabel("y [mm]")
    layer_text = axis.set_title(f"{title}\nlayer 01/57", fontsize=11)
    frames: list[Image.Image] = []
    for z in range(smoothed_crop.shape[0]):
        image.set_data(smoothed_crop[z])
        layer_text.set_text(f"{title}\nlayer {z + 1:02d}/57")
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=100)
        buffer.seek(0)
        frames.append(Image.open(buffer).convert("P", palette=Image.Palette.ADAPTIVE, colors=256))
        buffer.close()
    plt.close(figure)
    return frames


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--volume-glob", required=True)
    parser.add_argument("--prediction-glob", required=True)
    parser.add_argument(
        "--exclude-overlap-prediction-glob",
        help=(
            "Prediction HDF5 glob for a reference model. Candidates whose score-grid "
            "component overlaps a reference component are omitted."
        ),
    )
    parser.add_argument(
        "--only-overlap-prediction-glob",
        help=(
            "Prediction HDF5 glob for a reference model. Keep only candidates whose "
            "score-grid component overlaps a reference component."
        ),
    )
    parser.add_argument(
        "--exclude-overlap-threshold", type=float, default=0.90,
        help="Score threshold for target/reference component overlap (default: %(default).2f)",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--frame-duration-ms", type=int, default=160)
    parser.add_argument(
        "--display-crop-size",
        type=int,
        help="Displayed square crop in bins; default is the CNN crop size",
    )
    parser.add_argument(
        "--exclude-cell-id", type=int, action="append", default=[],
        help="Exclude a data/background cell; repeat for multiple cells",
    )
    args = parser.parse_args()
    if args.frame_duration_ms <= 0:
        raise ValueError("frame duration must be positive")
    if not 0.0 <= args.exclude_overlap_threshold <= 1.0:
        raise ValueError("exclude-overlap-threshold must be in [0, 1]")
    if (
        args.exclude_overlap_prediction_glob
        and args.only_overlap_prediction_glob
    ):
        parser.error("choose only one overlap selection mode")

    volume_paths = {path.name: path for path in map(Path, glob.glob(args.volume_glob))}
    prediction_paths = sorted(map(Path, glob.glob(args.prediction_glob)))
    if not volume_paths or not prediction_paths:
        raise ValueError("volume or prediction glob is empty")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    candidates: list[dict[str, object]] = []
    overlap_reference_glob = (
        args.exclude_overlap_prediction_glob or args.only_overlap_prediction_glob
    )
    reference_components = (
        components_by_event(overlap_reference_glob, args.exclude_overlap_threshold)
        if overlap_reference_glob else None
    )
    overlap_excluded = 0

    for prediction_path in prediction_paths:
        if prediction_path.name.endswith((".predictions.h5")):
            suffix = ".predictions.h5"
            volume_stem = prediction_path.name.removesuffix(suffix)
            volume_name = volume_stem + ".h5"
        elif prediction_path.name.endswith(("_t090.h5", "_t098.h5")):
            volume_stem = prediction_path.name.rsplit("_t", 1)[0]
            volume_name = volume_stem + ".h5"
        else:
            raise ValueError(
                f"prediction file {prediction_path.name} must end in "
                "'_predictions.h5' or '_t098.h5'"
            )
        if volume_name in volume_paths:
            volume_path = volume_paths[volume_name]
        else:
            # Predictions may omit the cell interval, while volume shards
            # retain it: data_scan_b000122_0000.h5 -> *_cellAAA-BBB.h5.
            matches = [p for n, p in volume_paths.items()
                       if n.startswith(volume_stem + "_cell")]
            if len(matches) != 1:
                raise ValueError(f"cannot uniquely match {prediction_path.name} to a source volume shard")
            volume_path = matches[0]
        with h5py.File(volume_path, "r") as source, h5py.File(prediction_path, "r") as prediction:
            group = prediction["candidates"]
            score_grids = prediction["window_presence_score"]
            x_positions = prediction["grid_x_start_bin"][:].astype(int)
            y_positions = prediction["grid_y_start_bin"][:].astype(int)
            crop_size = int(prediction.attrs["crop_size"])
            for index in range(len(group["event_id"])):
                event_index = int(group["event_index"][index])
                cell_id = int(source["cell_id"][event_index])
                if cell_id in set(args.exclude_cell_id):
                    continue
                grid_row = int(group["grid_row"][index])
                grid_col = int(group["grid_col"][index])
                if reference_components is not None:
                    component = component_at_score_threshold(
                        score_grids[event_index], args.exclude_overlap_threshold,
                        grid_row, grid_col,
                    )
                    references = reference_components.get(int(group["event_id"][index]), [])
                    overlaps_reference = any(component & reference for reference in references)
                    if args.exclude_overlap_prediction_glob and overlaps_reference:
                        overlap_excluded += 1
                        continue
                    if args.only_overlap_prediction_glob and not overlaps_reference:
                        overlap_excluded += 1
                        continue
                candidates.append({
                    "volume_path": volume_path,
                    "event_index": event_index,
                    "event_id": int(group["event_id"][index]),
                    "run_wall_brick": (
                        int(source["run_wall_brick"][event_index])
                        if "run_wall_brick" in source else -1
                    ),
                    "brick_name": (
                        str(source["brick_name"].asstr()[event_index])
                        if "brick_name" in source else ""
                    ),
                    "cell_id": cell_id,
                    "cell_x": int(source["cell_x"][event_index]),
                    "cell_y": int(source["cell_y"][event_index]),
                    "candidate_id": int(group["candidate_id"][index]),
                    "grid_row": grid_row,
                    "grid_col": grid_col,
                    "x_start": int(x_positions[grid_col]),
                    "y_start": int(y_positions[grid_row]),
                    "crop_size": crop_size,
                    "presence_score": float(group["presence_score"][index]),
                    "slope_x": float(group["slope_x"][index]),
                    "slope_y": float(group["slope_y"][index]),
                    "theta_pred_mrad": float(group["theta_pred_mrad"][index]),
                    "candidate_center_x_um": float(group["x_center_um"][index]),
                    "candidate_center_y_um": float(group["y_center_um"][index]),
                })

    candidates.sort(key=lambda row: float(row["presence_score"]), reverse=True)
    manifest_rows: list[dict[str, object]] = []
    for rank, candidate in enumerate(candidates, start=1):
        volume_path = Path(candidate["volume_path"])
        event_index = int(candidate["event_index"])
        x_start = int(candidate["x_start"])
        y_start = int(candidate["y_start"])
        crop_size = int(candidate["crop_size"])
        with h5py.File(volume_path, "r") as source:
            raw_volume = np.asarray(source["volumes_raw"][event_index])
            background_mu = float(source["background_mu"][event_index])
            x_edges = np.asarray(source["x_edges_um"][event_index])
            y_edges = np.asarray(source["y_edges_um"][event_index])
        display_crop_size = (
            crop_size if args.display_crop_size is None else args.display_crop_size
        )
        display_crop_size = fitted_display_crop_size(
            display_crop_size,
            crop_size,
            raw_volume.shape[1],
            raw_volume.shape[2],
        )
        display_x_start = centered_window_start(
            x_start, crop_size, display_crop_size, raw_volume.shape[2]
        )
        display_y_start = centered_window_start(
            y_start, crop_size, display_crop_size, raw_volume.shape[1]
        )
        # Smooth the complete 200x200 layer first, then center the requested
        # display window on the representative maximum-score CNN crop.
        smoothed = root_th2_smooth(raw_volume)
        crop = smoothed[
            :, display_y_start:display_y_start + display_crop_size,
            display_x_start:display_x_start + display_crop_size,
        ]
        raw_display_crop = raw_volume[
            :, display_y_start:display_y_start + display_crop_size,
            display_x_start:display_x_start + display_crop_size,
        ]
        x_centers_um = 0.5 * (
            x_edges[display_x_start:display_x_start + display_crop_size]
            + x_edges[display_x_start + 1:display_x_start + display_crop_size + 1]
        )
        y_centers_um = 0.5 * (
            y_edges[display_y_start:display_y_start + display_crop_size]
            + y_edges[display_y_start + 1:display_y_start + display_crop_size + 1]
        )
        start = estimate_shower_start(
            raw_display_crop,
            background_mu,
            x_centers_um,
            y_centers_um,
            predicted_slope_xy=np.asarray([
                candidate["slope_x"], candidate["slope_y"]
            ]),
            candidate_center_xy_um=np.asarray([
                candidate["candidate_center_x_um"], candidate["candidate_center_y_um"]
            ]),
        )
        extent_mm = (
            float(x_edges[display_x_start] / 1000.0),
            float(x_edges[display_x_start + display_crop_size] / 1000.0),
            float(y_edges[display_y_start] / 1000.0),
            float(y_edges[display_y_start + display_crop_size] / 1000.0),
        )
        tx_pred_mrad = 1000.0 * float(candidate["slope_x"]) / 27.0
        ty_pred_mrad = 1000.0 * float(candidate["slope_y"]) / 27.0
        brick_title = (
            f"{candidate['brick_name']} | " if str(candidate["brick_name"]) else ""
        )
        if start["start_layer"] is None:
            start_title = "start non trovato"
        else:
            start_symbol = "<=" if bool(start["start_at_first_observed_layer"]) else "="
            start_title = (
                f"start layer{start_symbol}{int(start['start_layer']):02d}/57 | "
                f"ROOT z-bin {int(start['start_root_z_bin']):02d} | "
                f"x={float(start['start_x_um']) / 1000.0:.3f} mm, "
                f"y={float(start['start_y_um']) / 1000.0:.3f} mm | "
                f"{start['start_quality']}"
            )
        title = (
            f"rank {rank:02d} | {brick_title}cell {int(candidate['cell_id'])} "
            f"(x={int(candidate['cell_x'])}, y={int(candidate['cell_y'])}) | "
            f"candidate {int(candidate['candidate_id'])}\n"
            f"P(signal)={float(candidate['presence_score']):.6f} | "
            f"theta_pred={float(candidate['theta_pred_mrad']):.2f} mrad\n"
            f"tx_pred={tx_pred_mrad:+.2f} mrad | ty_pred={ty_pred_mrad:+.2f} mrad\n"
            f"{start_title}"
        )
        frames = render_frames(crop, extent_mm=extent_mm, title=title)
        score_text = f"{float(candidate['presence_score']):.6f}".replace(".", "p")
        brick_filename = (
            f"{candidate['brick_name']}_" if str(candidate["brick_name"]) else ""
        )
        filename = (
            f"rank_{rank:03d}_{brick_filename}cell{int(candidate['cell_id']):03d}_"
            f"x{int(candidate['cell_x']):02d}_y{int(candidate['cell_y']):02d}_"
            f"cand{int(candidate['candidate_id']):02d}_score{score_text}.gif"
        )
        output = args.output_dir / filename
        durations = [1000] + [args.frame_duration_ms] * (len(frames) - 2) + [1000]
        frames[0].save(
            output,
            save_all=True,
            append_images=frames[1:],
            duration=durations,
            loop=0,
            optimize=False,
            disposal=2,
        )
        print(
            json.dumps({"status": "written", "rank": rank, "gif": str(output)}),
            flush=True,
        )
        manifest_rows.append({
            key: str(value) if isinstance(value, Path) else value
            for key, value in candidate.items()
        } | {
            "rank": rank,
            "gif": str(output),
            "extent_mm": list(extent_mm),
            "frames": len(frames),
            "display_crop_size": display_crop_size,
            "display_x_start": display_x_start,
            "display_y_start": display_y_start,
            "tx_pred_mrad": tx_pred_mrad,
            "ty_pred_mrad": ty_pred_mrad,
            **start,
        })

    manifest = {
        "candidate_count": len(candidates),
        "overlap_excluded_candidate_count": overlap_excluded,
        "overlap_reference_prediction_glob": overlap_reference_glob,
        "overlap_selection": (
            "exclude" if args.exclude_overlap_prediction_glob else
            "only" if args.only_overlap_prediction_glob else None
        ),
        "overlap_threshold": (
            args.exclude_overlap_threshold if overlap_reference_glob else None
        ),
        "excluded_cell_ids": sorted(set(args.exclude_cell_id)),
        "smoothing": "one ROOT TH2::Smooth(k5a)-equivalent pass on each complete 200x200 layer",
        "crop": (
            "display window centered on the representative 20x20 CNN window "
            "with maximum score in each connected cluster; clamped at cell boundaries"
        ),
        "display_crop_size": args.display_crop_size,
        "tagging_overlays": False,
        "frame_duration_ms": args.frame_duration_ms,
        "candidates": manifest_rows,
    }
    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    csv_path = args.output_dir / "candidate_starts.csv"
    csv_fields = (
        "rank", "brick_name", "event_id", "cell_id", "cell_x", "cell_y", "candidate_id",
        "presence_score", "theta_pred_mrad", "tx_pred_mrad", "ty_pred_mrad",
        "start_layer", "start_root_z_bin", "start_x_um", "start_y_um",
        "start_quality", "start_component_voxels", "start_component_layers",
        "start_component_charge", "component_selection",
        "selected_component_slope_x", "selected_component_slope_y",
        "poisson_hot_count_threshold", "gif",
    )
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(json.dumps({
        "candidate_count": len(candidates),
        "output_dir": str(args.output_dir.resolve()),
        "manifest": str(manifest_path.resolve()),
        "candidate_starts_csv": str(csv_path.resolve()),
    }, indent=2))


if __name__ == "__main__":
    main()
