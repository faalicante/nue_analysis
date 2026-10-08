#!/usr/bin/env python3
"""Export XZ/YZ q-hot projections with CNN and geometric slopes."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

from visualization.make_background_candidate_gifs import (
    centered_window_start,
    estimate_shower_start,
    fitted_display_crop_size,
    poisson_hot_threshold,
    root_th2_smooth,
)
from scanning.scan_cnn21d_volumes import cluster_contains_truth, connected_components


N_LAYERS = 57
PLATE_INDEX_OFFSET = 3
Z_STEP_UM = 1350.0


def visual_qhot_projections(
    raw_volume: np.ndarray,
    background_mu: float,
    x_start: int,
    y_start: int,
    size: int,
    *,
    smooth: bool,
    poisson_tail_probability: float,
) -> tuple[tuple[np.ndarray, np.ndarray], int]:
    """Build display-only XZ/YZ q-hot projections from raw or smoothed counts."""
    values = root_th2_smooth(raw_volume) if smooth else np.asarray(raw_volume, dtype=float)
    crop = values[:, y_start:y_start + size, x_start:x_start + size]
    threshold = poisson_hot_threshold(background_mu, poisson_tail_probability)
    significant = np.where(crop >= threshold, np.maximum(crop - background_mu, 0.0), 0.0)
    return (significant.sum(axis=1), significant.sum(axis=2)), threshold


def candidates_from_manifest(path: Path) -> list[dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return list(payload["candidates"])


def candidates_from_scan(
    volume_glob: str, prediction_glob: str, display_crop_size: int
) -> list[dict[str, object]]:
    volumes = {path.name: path for path in map(Path, glob.glob(volume_glob))}
    candidates: list[dict[str, object]] = []
    for prediction_path in sorted(map(Path, glob.glob(prediction_glob))):
        stem = prediction_path.stem
        for suffix in ("_t090", "_t098", "_predictions", ".predictions"):
            if stem.endswith(suffix):
                stem = stem.removesuffix(suffix)
                break
        volume_name = stem + ".h5"
        if volume_name in volumes:
            volume_path = volumes[volume_name]
        else:
            matches = [p for n, p in volumes.items() if n.startswith(stem + "_cell")]
            if len(matches) != 1:
                raise ValueError(f"cannot uniquely match {prediction_path.name} to a source volume")
            volume_path = matches[0]
        with h5py.File(prediction_path, "r") as prediction, h5py.File(volume_path, "r") as source:
            group = prediction["candidates"]
            grid_x = np.asarray(prediction["grid_x_start_bin"], dtype=int)
            grid_y = np.asarray(prediction["grid_y_start_bin"], dtype=int)
            crop_size = int(prediction.attrs["crop_size"])
            full_height = int(source["volumes_raw"].shape[-2])
            full_width = int(source["volumes_raw"].shape[-1])
            for index in range(len(group["event_id"])):
                event_index = int(group["event_index"][index])
                grid_col = int(group["grid_col"][index])
                grid_row = int(group["grid_row"][index])
                x_start = int(grid_x[grid_col])
                y_start = int(grid_y[grid_row])
                fitted_display_size = fitted_display_crop_size(
                    display_crop_size, crop_size, full_height, full_width
                )
                display_x = centered_window_start(
                    x_start, crop_size, fitted_display_size, full_width
                )
                display_y = centered_window_start(
                    y_start, crop_size, fitted_display_size, full_height
                )
                raw = np.asarray(source["volumes_raw"][
                    event_index, :, display_y:display_y + fitted_display_size,
                    display_x:display_x + fitted_display_size,
                ])
                x_edges = np.asarray(source["x_edges_um"][event_index])
                y_edges = np.asarray(source["y_edges_um"][event_index])
                x_centers = 0.5 * (
                    x_edges[display_x:display_x + fitted_display_size]
                    + x_edges[display_x + 1:display_x + fitted_display_size + 1]
                )
                y_centers = 0.5 * (
                    y_edges[display_y:display_y + fitted_display_size]
                    + y_edges[display_y + 1:display_y + fitted_display_size + 1]
                )
                slope_x = float(group["slope_x"][index])
                slope_y = float(group["slope_y"][index])
                estimate = estimate_shower_start(
                    raw,
                    float(source["background_mu"][event_index]),
                    x_centers,
                    y_centers,
                    predicted_slope_xy=np.asarray([slope_x, slope_y]),
                    candidate_center_xy_um=np.asarray([
                        group["x_center_um"][index], group["y_center_um"][index]
                    ]),
                )
                candidates.append({
                    "volume_path": str(volume_path),
                    "event_index": event_index,
                    "event_id": int(group["event_id"][index]),
                    "cell_id": int(source["cell_id"][event_index]),
                    "cell_x": int(source["cell_x"][event_index]),
                    "cell_y": int(source["cell_y"][event_index]),
                    "candidate_id": int(group["candidate_id"][index]),
                    "presence_score": float(group["presence_score"][index]),
                    "theta_pred_mrad": float(group["theta_pred_mrad"][index]),
                    "slope_x": slope_x,
                    "slope_y": slope_y,
                    "tx_pred_mrad": 1000.0 * slope_x / 27.0,
                    "ty_pred_mrad": 1000.0 * slope_y / 27.0,
                    "display_crop_size": fitted_display_size,
                    "display_x_start": display_x,
                    "display_y_start": display_y,
                    **estimate,
                })
    candidates.sort(key=lambda row: float(row["presence_score"]), reverse=True)
    for rank, candidate in enumerate(candidates, start=1):
        candidate["rank"] = rank
    return candidates


def event_lookup_from_volumes(volume_glob: str) -> dict[int, tuple[Path, int]]:
    lookup: dict[int, tuple[Path, int]] = {}
    for volume_path in sorted(map(Path, glob.glob(volume_glob))):
        with h5py.File(volume_path, "r") as source:
            for event_index, event_id in enumerate(np.asarray(source["event_id"], dtype=int)):
                if int(event_id) in lookup:
                    raise ValueError(f"event_id {int(event_id)} appears in multiple input volumes")
                lookup[int(event_id)] = (volume_path, int(event_index))
    if not lookup:
        raise ValueError(f"no input volumes matched {volume_glob!r}")
    return lookup


def truth_center_bins(
    source: h5py.File,
    event_index: int,
    truth_propagation_limit: int | None,
) -> tuple[float, float, float, float, dict[str, object]]:
    if truth_propagation_limit is None:
        center_x_um = float(source["truth/propagated_center_x_um"][event_index])
        center_y_um = float(source["truth/propagated_center_y_um"][event_index])
        center_x_bin = float(source["truth/propagated_center_x_bin"][event_index])
        center_y_bin = float(source["truth/propagated_center_y_bin"][event_index])
        metadata = {
            "truth_matching": "stored",
            "propagated_plates": None,
        }
    else:
        if truth_propagation_limit <= 0:
            raise ValueError("--truth-propagation-limit must be positive")
        crop_plate = int(source["truth/plate"][event_index]) - PLATE_INDEX_OFFSET
        propagated = min(truth_propagation_limit, N_LAYERS - crop_plate)
        factor = Z_STEP_UM * 0.5 * propagated / 1000.0
        center_x_um = (
            float(source["truth/xpos_um"][event_index])
            + float(source["truth/ltx_mrad_corrected"][event_index]) * factor
        )
        center_y_um = (
            float(source["truth/ypos_um"][event_index])
            + float(source["truth/lty_mrad_corrected"][event_index]) * factor
        )
        x_centers = 0.5 * (
            np.asarray(source["x_edges_um"][event_index, :-1], dtype=float)
            + np.asarray(source["x_edges_um"][event_index, 1:], dtype=float)
        )
        y_centers = 0.5 * (
            np.asarray(source["y_edges_um"][event_index, :-1], dtype=float)
            + np.asarray(source["y_edges_um"][event_index, 1:], dtype=float)
        )
        center_x_bin = float(np.interp(center_x_um, x_centers, np.arange(x_centers.size)))
        center_y_bin = float(np.interp(center_y_um, y_centers, np.arange(y_centers.size)))
        metadata = {
            "truth_matching": f"recomputed_propagation_limit_{truth_propagation_limit}",
            "crop_plate": crop_plate,
            "propagated_plates": propagated,
        }
    theta_true = float(1000.0 * math.atan(
        math.hypot(
            float(source["truth/ltx_mrad_corrected"][event_index]),
            float(source["truth/lty_mrad_corrected"][event_index]),
        ) / 1000.0
    ))
    metadata.update({
        "truth_center_x_um": center_x_um,
        "truth_center_y_um": center_y_um,
        "truth_vertex_x_um": float(source["truth/xpos_um"][event_index]),
        "truth_vertex_y_um": float(source["truth/ypos_um"][event_index]),
        "truth_center_x_bin": center_x_bin,
        "truth_center_y_bin": center_y_bin,
        "truth_plate": int(source["truth/plate"][event_index]),
        "truth_vertex_layer": int(source["truth/plate"][event_index]) - PLATE_INDEX_OFFSET,
        "theta_true_mrad": theta_true,
        "tx_true_mrad": float(source["truth/ltx_mrad_corrected"][event_index]),
        "ty_true_mrad": float(source["truth/lty_mrad_corrected"][event_index]),
    })
    return center_x_bin, center_y_bin, center_x_um, center_y_um, metadata


def nearest_truth_window(
    truth_x_bin: float,
    truth_y_bin: float,
    x_positions: np.ndarray,
    y_positions: np.ndarray,
    crop_size: int,
) -> tuple[int, int, float]:
    x_centers = x_positions.astype(float) + crop_size / 2.0
    y_centers = y_positions.astype(float) + crop_size / 2.0
    distances = np.hypot(
        y_centers[:, None] - truth_y_bin,
        x_centers[None, :] - truth_x_bin,
    )
    row, col = np.unravel_index(int(np.argmin(distances)), distances.shape)
    return int(row), int(col), float(distances[row, col])


def nearest_truth_location(
    locations: list[tuple[int, int]],
    truth_x_bin: float,
    truth_y_bin: float,
    x_positions: np.ndarray,
    y_positions: np.ndarray,
    crop_size: int,
) -> tuple[int, int, float]:
    if not locations:
        raise ValueError("cannot choose a nearest location from an empty list")
    best_row, best_col = min(
        locations,
        key=lambda location: math.hypot(
            y_positions[location[0]] + crop_size / 2.0 - truth_y_bin,
            x_positions[location[1]] + crop_size / 2.0 - truth_x_bin,
        ),
    )
    distance = math.hypot(
        y_positions[best_row] + crop_size / 2.0 - truth_y_bin,
        x_positions[best_col] + crop_size / 2.0 - truth_x_bin,
    )
    return int(best_row), int(best_col), float(distance)


def window_contains_truth(
    row: int,
    col: int,
    truth_x_bin: float,
    truth_y_bin: float,
    x_positions: np.ndarray,
    y_positions: np.ndarray,
    crop_size: int,
) -> bool:
    return (
        x_positions[col] <= truth_x_bin < x_positions[col] + crop_size
        and y_positions[row] <= truth_y_bin < y_positions[row] + crop_size
    )


def signal_candidates_from_scan(
    volume_glob: str,
    prediction_path: Path,
    display_crop_size: int,
    *,
    threshold: float | None = None,
    truth_propagation_limit: int | None = None,
    selection: str = "all",
    theta_min: float | None = None,
    found_window: str = "truth-containing-nearest",
) -> list[dict[str, object]]:
    event_lookup = event_lookup_from_volumes(volume_glob)
    candidates: list[dict[str, object]] = []
    with h5py.File(prediction_path, "r") as prediction:
        scores = np.asarray(prediction["window_presence_score"], dtype=float)
        slopes = np.asarray(prediction["window_slope_xy"], dtype=float)
        prediction_event_ids = np.asarray(prediction["event_id"], dtype=int)
        grid_x = np.asarray(prediction["grid_x_start_bin"], dtype=int)
        grid_y = np.asarray(prediction["grid_y_start_bin"], dtype=int)
        crop_size = int(prediction.attrs["crop_size"])
        threshold_value = float(
            prediction.attrs["presence_threshold"] if threshold is None else threshold
        )
        for prediction_event_index, event_id in enumerate(prediction_event_ids):
            if int(event_id) not in event_lookup:
                raise ValueError(f"event_id {int(event_id)} not found in {volume_glob!r}")
            volume_path, source_event_index = event_lookup[int(event_id)]
            with h5py.File(volume_path, "r") as source:
                full_size = int(source["volumes_raw"].shape[-1])
                truth_x_bin, truth_y_bin, truth_x_um, truth_y_um, truth = truth_center_bins(
                    source, source_event_index, truth_propagation_limit
                )
                if theta_min is not None and float(truth["theta_true_mrad"]) < theta_min:
                    continue
                components = connected_components(scores[prediction_event_index] >= threshold_value)
                matched_components = [
                    component for component in components
                    if cluster_contains_truth(
                        component, truth_x_bin, truth_y_bin,
                        grid_x.tolist(), grid_y.tolist(), crop_size,
                    )
                ]
                if matched_components:
                    if selection in ("all", "found"):
                        for candidate_id, component in enumerate(matched_components):
                            if found_window == "cluster-max-score":
                                row, col = max(
                                    component,
                                    key=lambda location: float(scores[prediction_event_index][location]),
                                )
                                distance = 0.0
                            else:
                                truth_containing = [
                                    location for location in component
                                    if window_contains_truth(
                                        location[0], location[1],
                                        truth_x_bin, truth_y_bin,
                                        grid_x, grid_y, crop_size,
                                    )
                                ]
                                row, col, distance = nearest_truth_location(
                                    truth_containing or component,
                                    truth_x_bin, truth_y_bin,
                                    grid_x, grid_y, crop_size,
                                )
                            candidates.append(signal_candidate_row(
                                source, source_event_index, volume_path,
                                grid_x, grid_y, crop_size, display_crop_size, full_size,
                                row, col,
                                presence_score=float(scores[prediction_event_index, row, col]),
                                slope=slopes[prediction_event_index, row, col],
                                candidate_id=candidate_id,
                                kind="found",
                                threshold=threshold_value,
                                truth=truth,
                                nearest_distance_bins=distance,
                                found_window=found_window,
                                truth_x_um=truth_x_um,
                                truth_y_um=truth_y_um,
                            ))
                elif selection in ("all", "missed"):
                    row, col, distance = nearest_truth_window(
                        truth_x_bin, truth_y_bin, grid_x, grid_y, crop_size
                    )
                    candidates.append(signal_candidate_row(
                        source, source_event_index, volume_path,
                        grid_x, grid_y, crop_size, display_crop_size, full_size,
                        row, col,
                        presence_score=float(scores[prediction_event_index, row, col]),
                        slope=slopes[prediction_event_index, row, col],
                        candidate_id=0,
                        kind="missed",
                        threshold=threshold_value,
                        truth=truth,
                        nearest_distance_bins=distance,
                        found_window=found_window,
                        truth_x_um=truth_x_um,
                        truth_y_um=truth_y_um,
                    ))
    candidates.sort(
        key=lambda row: (
            0 if row["signal_scan_kind"] == "found" else 1,
            -float(row.get("theta_true_mrad", 0.0)),
            -float(row["presence_score"]),
        )
    )
    kind_counts = {"found": 0, "missed": 0}
    for candidate in candidates:
        kind = str(candidate["signal_scan_kind"])
        kind_counts[kind] += 1
        candidate["rank"] = kind_counts[kind]
    return candidates


def signal_candidate_row(
    source: h5py.File,
    event_index: int,
    volume_path: Path,
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    crop_size: int,
    display_crop_size: int,
    full_size: int,
    row: int,
    col: int,
    *,
    presence_score: float,
    slope: np.ndarray,
    candidate_id: int,
    kind: str,
    threshold: float,
    truth: dict[str, object],
    nearest_distance_bins: float,
    found_window: str,
    truth_x_um: float,
    truth_y_um: float,
) -> dict[str, object]:
    x_start = int(grid_x[col])
    y_start = int(grid_y[row])
    display_x = centered_window_start(x_start, crop_size, display_crop_size, full_size)
    display_y = centered_window_start(y_start, crop_size, display_crop_size, full_size)
    raw = np.asarray(source["volumes_raw"][
        event_index, :, display_y:display_y + display_crop_size,
        display_x:display_x + display_crop_size,
    ])
    x_edges = np.asarray(source["x_edges_um"][event_index])
    y_edges = np.asarray(source["y_edges_um"][event_index])
    x_centers = 0.5 * (
        x_edges[display_x:display_x + display_crop_size]
        + x_edges[display_x + 1:display_x + display_crop_size + 1]
    )
    y_centers = 0.5 * (
        y_edges[display_y:display_y + display_crop_size]
        + y_edges[display_y + 1:display_y + display_crop_size + 1]
    )
    slope_x = float(slope[0])
    slope_y = float(slope[1])
    estimate = estimate_shower_start(
        raw,
        float(source["background_mu"][event_index]),
        x_centers,
        y_centers,
        predicted_slope_xy=np.asarray([slope_x, slope_y]),
        candidate_center_xy_um=np.asarray([truth_x_um, truth_y_um]),
    )
    return {
        "volume_path": str(volume_path),
        "event_index": event_index,
        "event_id": int(source["event_id"][event_index]),
        "cell_id": int(source["cell_id"][event_index]) if "cell_id" in source else -1,
        "cell_x": int(source["cell_x"][event_index]) if "cell_x" in source else -1,
        "cell_y": int(source["cell_y"][event_index]) if "cell_y" in source else -1,
        "candidate_id": candidate_id,
        "presence_score": presence_score,
        "theta_pred_mrad": float(1000.0 * math.atan(math.hypot(slope_x, slope_y) / 27.0)),
        "slope_x": slope_x,
        "slope_y": slope_y,
        "tx_pred_mrad": 1000.0 * slope_x / 27.0,
        "ty_pred_mrad": 1000.0 * slope_y / 27.0,
        "display_crop_size": display_crop_size,
        "display_x_start": display_x,
        "display_y_start": display_y,
        "scan_grid_row": int(row),
        "scan_grid_col": int(col),
        "scan_crop_x_start": x_start,
        "scan_crop_y_start": y_start,
        "scan_threshold": threshold,
        "signal_scan_kind": kind,
        "found_window": found_window,
        "nearest_truth_window_distance_bins": nearest_distance_bins,
        **truth,
        **estimate,
    }


def slope_line(start_um: float, slope_mrad: float, start_layer: int) -> np.ndarray:
    layers = np.arange(1, 58, dtype=float)
    return start_um + slope_mrad * Z_STEP_UM * (layers - start_layer) / 1000.0


def render_candidate(
    candidate: dict[str, object],
    output: Path,
    sample_name: str,
    *,
    smooth_projections: bool,
    projection_poisson_tail_probability: float,
) -> dict[str, object]:
    with h5py.File(str(candidate["volume_path"]), "r") as source:
        event_index = int(candidate["event_index"])
        x0 = int(candidate["display_x_start"])
        y0 = int(candidate["display_y_start"])
        size = int(candidate["display_crop_size"])
        raw_volume = np.asarray(source["volumes_raw"][event_index], dtype=float)
        mu = float(source["background_mu"][event_index])
        x_edges = np.asarray(source["x_edges_um"][event_index, x0:x0 + size + 1], dtype=float)
        y_edges = np.asarray(source["y_edges_um"][event_index, y0:y0 + size + 1], dtype=float)

    projections, projection_threshold = visual_qhot_projections(
        raw_volume,
        mu,
        x0,
        y0,
        size,
        smooth=smooth_projections,
        poisson_tail_probability=projection_poisson_tail_probability,
    )
    physical_edges = (x_edges / 1000.0, y_edges / 1000.0)
    bin_widths = (float(np.median(np.diff(x_edges))), float(np.median(np.diff(y_edges))))
    cnn_slopes = (float(candidate["tx_pred_mrad"]), float(candidate["ty_pred_mrad"]))
    component_slope_x = candidate.get("selected_component_slope_x")
    component_slope_y = candidate.get("selected_component_slope_y")
    component_slopes = (
        (
            1000.0 * float(component_slope_x) * bin_widths[0] / Z_STEP_UM
            if component_slope_x is not None else float("nan")
        ),
        (
            1000.0 * float(component_slope_y) * bin_widths[1] / Z_STEP_UM
            if component_slope_y is not None else float("nan")
        ),
    )
    fallback_start = (
        float(candidate.get("truth_center_x_um", np.nan)),
        float(candidate.get("truth_center_y_um", np.nan)),
    )
    fallback_layer = int(candidate.get("truth_plate", 1))
    start = (
        float(candidate["start_x_um"]) if candidate.get("start_x_um") is not None else fallback_start[0],
        float(candidate["start_y_um"]) if candidate.get("start_y_um") is not None else fallback_start[1],
    )
    start_layer = (
        int(candidate["start_layer"]) if candidate.get("start_layer") is not None else fallback_layer
    )
    component_theta = float(
        1000.0 * math.atan(math.hypot(*component_slopes) / 1000.0)
    )
    true_slopes = (
        candidate.get("tx_true_mrad"),
        candidate.get("ty_true_mrad"),
    )
    truth_vertex = (
        candidate.get("truth_vertex_x_um"),
        candidate.get("truth_vertex_y_um"),
        candidate.get("truth_vertex_layer"),
    )
    has_truth_vertex = all(value is not None for value in truth_vertex)

    figure, axes = plt.subplots(1, 2, figsize=(12.4, 5.2), constrained_layout=True)
    for axis, projection, edges, coordinate, cnn_slope, fit_slope, true_slope, start_coordinate in zip(
        axes, projections, physical_edges, ("x", "y"), cnn_slopes,
        component_slopes, true_slopes, start, strict=True,
    ):
        positive = projection[projection > 0]
        vmax = float(np.quantile(positive, 0.995)) if positive.size else 1.0
        # raw is indexed [z, y, x].  Plot with z on the horizontal axis
        # (XZ/YZ, i.e. swap the previous coordinate-vs-layer orientation).
        mesh = axis.pcolormesh(
            np.arange(0.5, 58.5), edges, projection.T, cmap="viridis",
            norm=PowerNorm(gamma=0.55, vmin=0.0, vmax=max(vmax, 1.0)), shading="flat",
        )
        layers = np.arange(1, 58)
        axis.plot(
            layers, slope_line(start_coordinate, cnn_slope, start_layer) / 1000.0,
            color="#ff3b30", lw=2, label=f"CNN {cnn_slope:+.1f} mrad",
        )
        if np.isfinite(fit_slope):
            axis.plot(
                layers, slope_line(start_coordinate, fit_slope, start_layer) / 1000.0,
                color="white", lw=1.8, ls="--", label=f"fit {fit_slope:+.1f} mrad",
            )
        if true_slope is not None:
            axis.plot(
                layers,
                slope_line(
                    float(truth_vertex[0 if coordinate == "x" else 1]),
                    float(true_slope),
                    int(truth_vertex[2]),
                ) / 1000.0 if has_truth_vertex else
                slope_line(start_coordinate, float(true_slope), start_layer) / 1000.0,
                color="#2ed573", lw=1.7, ls=":", label=f"MC slope {float(true_slope):+.1f} mrad",
            )
        if has_truth_vertex:
            vertex_coordinate = float(truth_vertex[0 if coordinate == "x" else 1])
            axis.scatter(
                [int(truth_vertex[2])], [vertex_coordinate / 1000.0],
                color="#00e676", marker="*", s=105, edgecolors="black",
                linewidths=0.7, zorder=5, label="MC vertice originale",
            )
        axis.scatter(
            [start_layer], [start_coordinate / 1000.0], color="#ffcc00",
            marker="x", s=44, lw=2.2, label="inizio stimato",
        )
        axis.set_xlim(0.5, 57.5)
        axis.set_ylim(edges[0], edges[-1])
        axis.set_xlabel("layer XYPseg (z)")
        axis.set_ylabel(f"{coordinate} [mm]")
        axis.set_title(f"{coordinate.upper()}Z")
        axis.grid(alpha=0.16)
        axis.legend(loc="best", fontsize=8.5, framealpha=0.84)
        display_kind = "smooth q-hot" if smooth_projections else "q-hot"
        figure.colorbar(mesh, ax=axis, pad=0.01, label=f"eccesso {display_kind} integrato")
    kind = candidate.get("signal_scan_kind")
    kind_label = f"{kind} " if kind else ""
    theta_true = (
        f" · θ MC={float(candidate['theta_true_mrad']):.1f} mrad"
        if "theta_true_mrad" in candidate else ""
    )
    cell_label = f"cell {candidate['cell_id']} · " if int(candidate["cell_id"]) >= 0 else ""
    figure.suptitle(
        f"{sample_name} {kind_label}rank {candidate['rank']} · event {candidate.get('event_id', -1)} · "
        f"{cell_label}score={float(candidate['presence_score']):.3f} · "
        f"θ CNN={float(candidate['theta_pred_mrad']):.1f} mrad · θ fit={component_theta:.1f} mrad"
        f"{theta_true}",
        fontsize=13,
    )
    figure.savefig(output, dpi=135, bbox_inches="tight")
    plt.close(figure)
    row = {
        "rank": int(candidate["rank"]),
        "sample": sample_name,
        "signal_scan_kind": candidate.get("signal_scan_kind"),
        "event_id": int(candidate.get("event_id", -1)),
        "cell_id": int(candidate["cell_id"]),
        "cell_x": int(candidate["cell_x"]),
        "cell_y": int(candidate["cell_y"]),
        "candidate_id": int(candidate["candidate_id"]),
        "presence_score": float(candidate["presence_score"]),
        "theta_cnn_mrad": float(candidate["theta_pred_mrad"]),
        "tx_cnn_mrad": cnn_slopes[0],
        "ty_cnn_mrad": cnn_slopes[1],
        "theta_fit_mrad": component_theta,
        "tx_fit_mrad": component_slopes[0],
        "ty_fit_mrad": component_slopes[1],
        "slope_difference_mrad": (
            float(math.dist(cnn_slopes, component_slopes))
            if all(np.isfinite(component_slopes)) else None
        ),
        "start_layer": start_layer,
        "start_x_um": start[0],
        "start_y_um": start[1],
        "start_quality": candidate["start_quality"],
        "projection_smoothed": smooth_projections,
        "projection_poisson_tail_probability": projection_poisson_tail_probability,
        "projection_poisson_hot_count_threshold": projection_threshold,
        "image": str(output.resolve()),
    }
    for key in (
        "scan_threshold", "scan_grid_row", "scan_grid_col",
        "scan_crop_x_start", "scan_crop_y_start",
        "found_window",
        "nearest_truth_window_distance_bins",
        "truth_matching", "truth_plate", "crop_plate", "propagated_plates",
        "theta_true_mrad", "tx_true_mrad", "ty_true_mrad",
        "truth_center_x_um", "truth_center_y_um",
        "truth_vertex_x_um", "truth_vertex_y_um", "truth_vertex_layer",
        "truth_center_x_bin", "truth_center_y_bin",
    ):
        if key in candidate:
            row[key] = candidate[key]
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--manifest", type=Path)
    source.add_argument("--volume-glob")
    source.add_argument("--signal-volume-glob")
    parser.add_argument("--prediction-glob")
    parser.add_argument("--signal-predictions", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--display-crop-size", type=int, default=40)
    parser.add_argument(
        "--smooth-projections",
        action="store_true",
        help="Use one ROOT TH2::Smooth(k5a)-equivalent pass for the displayed q-hot maps only.",
    )
    parser.add_argument(
        "--projection-poisson-tail-probability",
        type=float,
        default=1.0e-4,
        help=(
            "Poisson upper-tail probability for displayed q-hot maps only "
            "(larger value means a lower count threshold; default: %(default)g)."
        ),
    )
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--truth-propagation-limit", type=int)
    parser.add_argument("--signal-selection", choices=("all", "found", "missed"), default="all")
    parser.add_argument(
        "--found-window",
        choices=("truth-containing-nearest", "cluster-max-score"),
        default="truth-containing-nearest",
        help=(
            "For found signal events, draw either the truth-containing matched "
            "window nearest to MC or the cluster's maximum-score window."
        ),
    )
    parser.add_argument("--theta-min", type=float)
    args = parser.parse_args()
    if not 0.0 < args.projection_poisson_tail_probability < 1.0:
        parser.error("--projection-poisson-tail-probability must be between 0 and 1")
    if args.manifest:
        candidates = candidates_from_manifest(args.manifest)
    elif args.volume_glob:
        if not args.prediction_glob:
            parser.error("--prediction-glob is required with --volume-glob")
        candidates = candidates_from_scan(
            args.volume_glob, args.prediction_glob, args.display_crop_size
        )
    else:
        if not args.signal_predictions:
            parser.error("--signal-predictions is required with --signal-volume-glob")
        candidates = signal_candidates_from_scan(
            args.signal_volume_glob,
            args.signal_predictions,
            args.display_crop_size,
            threshold=args.threshold,
            truth_propagation_limit=args.truth_propagation_limit,
            selection=args.signal_selection,
            theta_min=args.theta_min,
            found_window=args.found_window,
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        prefix = (
            f"{candidate['signal_scan_kind']}_" if candidate.get("signal_scan_kind") else ""
        )
        cell_part = (
            f"cell_{int(candidate['cell_id']):03d}_"
            if int(candidate["cell_id"]) >= 0 else ""
        )
        output = args.output_dir / (
            f"{prefix}rank_{int(candidate['rank']):03d}_event_{int(candidate.get('event_id', -1)):05d}_"
            f"{cell_part}"
            f"score_{float(candidate['presence_score']):.6f}_xz_yz.png"
        )
        rows.append(render_candidate(
            candidate,
            output,
            args.sample_name,
            smooth_projections=args.smooth_projections,
            projection_poisson_tail_probability=args.projection_poisson_tail_probability,
        ))
        print(json.dumps({"rank": candidate["rank"], "image": str(output)}), flush=True)
    manifest_csv = args.output_dir / "manifest.csv"
    if rows:
        fieldnames = sorted({key for row in rows for key in row})
        with manifest_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    (args.output_dir / "manifest.json").write_text(
        json.dumps({"count": len(rows), "sample": args.sample_name, "candidates": rows}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"count": len(rows), "manifest": str(manifest_csv.resolve())}, indent=2))


if __name__ == "__main__":
    main()
