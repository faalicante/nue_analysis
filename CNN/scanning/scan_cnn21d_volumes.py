#!/usr/bin/env python3
"""Scan 57x200x200 HDF5 volumes with overlapping CNN crops and cluster detections."""

from __future__ import annotations

import argparse
import json
from collections import deque
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import torch

from training.cnn_dataset import local_mpv_sigma_per_z, normalize_per_slice
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, choose_device, load_config


DEFAULT_THRESHOLD = 0.208049014210701
DEFAULT_TRUTH_PROPAGATION_LIMIT = 40
PLATE_INDEX_OFFSET = 3
Z_BIN_SIZE_UM = 1350.0


def sliding_positions(length: int, crop_size: int, stride: int) -> list[int]:
    if crop_size <= 0 or stride <= 0 or crop_size > length:
        raise ValueError("invalid length, crop size, or stride")
    positions = list(range(0, length - crop_size + 1, stride))
    if positions[-1] != length - crop_size:
        positions.append(length - crop_size)
    return positions


def connected_components(mask: np.ndarray) -> list[list[tuple[int, int]]]:
    """Eight-connected components in the crop-grid plane."""
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 2:
        raise ValueError("component mask must be two-dimensional")
    seen = np.zeros_like(mask)
    components: list[list[tuple[int, int]]] = []
    for row, col in zip(*np.nonzero(mask), strict=True):
        if seen[row, col]:
            continue
        queue = deque([(int(row), int(col))])
        seen[row, col] = True
        component: list[tuple[int, int]] = []
        while queue:
            current_row, current_col = queue.popleft()
            component.append((current_row, current_col))
            for delta_row in (-1, 0, 1):
                for delta_col in (-1, 0, 1):
                    if delta_row == 0 and delta_col == 0:
                        continue
                    next_row, next_col = current_row + delta_row, current_col + delta_col
                    if (
                        0 <= next_row < mask.shape[0] and 0 <= next_col < mask.shape[1]
                        and mask[next_row, next_col] and not seen[next_row, next_col]
                    ):
                        seen[next_row, next_col] = True
                        queue.append((next_row, next_col))
        components.append(component)
    return components


def crop_batches(
    volume: np.ndarray,
    background_mu: float,
    y_positions: list[int],
    x_positions: list[int],
    crop_size: int,
    batch_size: int,
    normalization_mode: str,
) -> Iterable[tuple[list[tuple[int, int]], torch.Tensor]]:
    locations = [(row, col) for row in range(len(y_positions)) for col in range(len(x_positions))]
    for first in range(0, len(locations), batch_size):
        batch_locations = locations[first : first + batch_size]
        crops = np.stack([
            normalize_scan_crop(
                volume, background_mu, y_positions[row], x_positions[col], crop_size,
                normalization_mode,
            )
            for row, col in batch_locations
        ])
        yield batch_locations, torch.from_numpy(crops[:, None]).float()


def normalize_scan_crop(
    volume: np.ndarray,
    background_mu: float,
    y_start: int,
    x_start: int,
    crop_size: int,
    normalization_mode: str,
) -> np.ndarray:
    """Normalize a scan window with the same local reference as training.

    The training crops carry a 32x32 source region around their 20x20 network
    crop.  During scanning we rebuild that 32x32 region around every window;
    at the 200x200 borders it is shifted inward without padding.
    """

    crop = volume[:, y_start : y_start + crop_size, x_start : x_start + crop_size]
    if normalization_mode.startswith("mpv_z"):
        reference_size = min(32, volume.shape[-2], volume.shape[-1])
        center_y = y_start + crop_size // 2
        center_x = x_start + crop_size // 2
        reference_y = min(max(center_y - reference_size // 2, 0), volume.shape[-2] - reference_size)
        reference_x = min(max(center_x - reference_size // 2, 0), volume.shape[-1] - reference_size)
        reference = volume[
            :, reference_y : reference_y + reference_size,
            reference_x : reference_x + reference_size,
        ]
        mpv_z, sigma_z = local_mpv_sigma_per_z(reference)
        return normalize_per_slice(
            crop, background_mu, mode=normalization_mode, mpv_z=mpv_z, sigma_z=sigma_z
        )
    return normalize_per_slice(crop, background_mu, mode=normalization_mode)


def augment_volume_background(
    volume: np.ndarray,
    background_mu: float,
    target_mu: float,
    *,
    seed: int,
) -> tuple[np.ndarray, float, float]:
    """Raise one scan volume to target mu while preserving raw zero voxels."""

    if not np.isfinite(background_mu) or background_mu <= 0.0:
        raise ValueError("background_mu must be finite and positive")
    if not np.isfinite(target_mu) or target_mu <= 0.0:
        raise ValueError("Poisson background target mu must be finite and positive")
    effective_mu = max(float(background_mu), float(target_mu))
    added_mu = effective_mu - float(background_mu)
    if added_mu <= 0.0:
        return volume, effective_mu, 0.0
    rng = np.random.default_rng(seed)
    noise = rng.poisson(added_mu, size=volume.shape).astype(np.float32)
    noise[np.asarray(volume) <= 0] = 0.0
    return np.asarray(volume, dtype=np.float32) + noise, effective_mu, added_mu


def local_crop_center(edges: np.ndarray, start: int, crop_size: int) -> float:
    return float(0.5 * (edges[start] + edges[start + crop_size]))


def coordinate_to_local_bin(edges: np.ndarray, value: float) -> float:
    centers = 0.5 * (edges[:-1] + edges[1:])
    return float(np.interp(value, centers, np.arange(centers.size, dtype=np.float64)))


def truth_reference_bins(
    source: h5py.File,
    event_index: int,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    plates: int,
    propagation_limit: int,
) -> tuple[float, float]:
    """Return the MC reference point used for truth matching.

    The scan volume itself is unchanged. This only selects where along the
    MC track the truth point is placed for the cluster-matching diagnostic.
    """
    if propagation_limit <= 0:
        raise ValueError("truth propagation limit must be positive")
    truth = source["truth"]
    crop_plate = int(truth["plate"][event_index]) - PLATE_INDEX_OFFSET
    propagated_plates = min(propagation_limit, plates - crop_plate)
    if propagated_plates < 0:
        raise ValueError(
            f"event {event_index}: crop plate {crop_plate} is beyond the scan volume"
        )
    half_length_um = 0.5 * propagated_plates * Z_BIN_SIZE_UM
    slope_factor = half_length_um / 1000.0
    x_um = float(truth["xpos_um"][event_index]) + float(
        truth["ltx_mrad_corrected"][event_index]
    ) * slope_factor
    y_um = float(truth["ypos_um"][event_index]) + float(
        truth["lty_mrad_corrected"][event_index]
    ) * slope_factor
    return coordinate_to_local_bin(x_edges, x_um), coordinate_to_local_bin(y_edges, y_um)


def cluster_contains_truth(
    component: list[tuple[int, int]],
    truth_x_bin: float,
    truth_y_bin: float,
    x_positions: list[int],
    y_positions: list[int],
    crop_size: int,
) -> bool:
    return any(
        x_positions[col] <= truth_x_bin < x_positions[col] + crop_size
        and y_positions[row] <= truth_y_bin < y_positions[row] + crop_size
        for row, col in component
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument(
        "--regression-checkpoint", type=Path,
        help="Optional separate checkpoint used only for slope_xy.",
    )
    parser.add_argument(
        "--regression-min-score",
        type=float,
        help=(
            "With a separate regression checkpoint, evaluate it only on windows at or "
            "above this classification score; other stored slopes are NaN."
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--crop-size", type=int, default=20)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument(
        "--device",
        default="auto",
        choices=("auto", "cpu", "cuda", "mps"),
        help="Inference device; auto selects CUDA, then MPS, then CPU.",
    )
    parser.add_argument(
        "--poisson-background-target-mu",
        type=float,
        help="Deterministically raise each input volume to this background mu before scanning.",
    )
    parser.add_argument(
        "--poisson-background-seed",
        type=int,
        default=20260910,
        help="Base seed for deterministic scan-time Poisson background augmentation.",
    )
    parser.add_argument(
        "--truth-propagation-limit",
        type=int,
        default=DEFAULT_TRUTH_PROPAGATION_LIMIT,
        help="Maximum number of z plates used for the MC truth-matching point.",
    )
    args = parser.parse_args()
    if args.poisson_background_seed < 0:
        raise ValueError("--poisson-background-seed must be non-negative")
    if args.regression_min_score is not None:
        if args.regression_checkpoint is None:
            raise ValueError("--regression-min-score requires --regression-checkpoint")
        if not 0.0 <= args.regression_min_score <= args.threshold:
            raise ValueError("--regression-min-score must be in [0, threshold]")

    config = load_config(args.config)
    normalization_mode = str(
        config.get("data", {}).get("dataset", {}).get(
            "normalization_mode", "background_mu_sigma"
        )
    )
    device = choose_device(args.device)
    model = build_model(config).to(device)
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    regression_model = model
    if args.regression_checkpoint is not None:
        regression_model = build_model(config).to(device)
        regression_checkpoint = torch.load(
            args.regression_checkpoint, map_location="cpu", weights_only=False
        )
        regression_model.load_state_dict(regression_checkpoint["model_state"])
        regression_model.eval()

    candidates: dict[str, list[object]] = {
        key: [] for key in (
            "event_index", "event_id", "candidate_id", "grid_row", "grid_col",
            "cluster_size", "presence_score", "slope_x", "slope_y",
            "theta_pred_mrad", "phi_pred_rad", "x_center_um", "y_center_um",
            "truth_matched",
        )
    }
    event_reports: list[dict[str, object]] = []
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.input, "r") as source, h5py.File(args.output, "w") as output:
        volumes = source["volumes_raw"]
        n, plates, height, width = volumes.shape
        if height < args.crop_size or width < args.crop_size:
            raise ValueError(
                f"input volume is {height}x{width} bins, smaller than the "
                f"requested {args.crop_size}x{args.crop_size} CNN crop"
            )
        x_positions = sliding_positions(width, args.crop_size, args.stride)
        y_positions = sliding_positions(height, args.crop_size, args.stride)
        scores_ds = output.create_dataset(
            "window_presence_score", shape=(n, len(y_positions), len(x_positions)), dtype=np.float32
        )
        slopes_ds = output.create_dataset(
            "window_slope_xy", shape=(n, len(y_positions), len(x_positions), 2), dtype=np.float32
        )
        output.create_dataset("event_id", data=source["event_id"][:])
        # Preserve optional background-cell provenance so candidate rows can be
        # mapped back to their input ROOT without decoding event_id.
        for metadata_name in (
            "cell_id", "cell_x", "cell_y", "source_file",
            "run_wall_brick", "rwb_token", "brick_name",
        ):
            if metadata_name in source:
                source.copy(metadata_name, output)
        output.create_dataset("grid_x_start_bin", data=np.asarray(x_positions, dtype=np.int16))
        output.create_dataset("grid_y_start_bin", data=np.asarray(y_positions, dtype=np.int16))
        output.attrs["crop_size"] = args.crop_size
        output.attrs["stride"] = args.stride
        output.attrs["presence_threshold"] = args.threshold
        output.attrs["truth_propagation_limit"] = args.truth_propagation_limit
        output.attrs["checkpoint"] = str(args.checkpoint)
        output.attrs["normalization_mode"] = normalization_mode
        if normalization_mode.startswith("mpv_z"):
            output.attrs["local_normalization_reference"] = "32x32 region centered on each 20x20 scan window"
        output.attrs["regression_checkpoint"] = str(
            args.regression_checkpoint or args.checkpoint
        )
        if args.regression_min_score is not None:
            output.attrs["regression_min_score"] = args.regression_min_score
        if args.poisson_background_target_mu is not None:
            output.attrs["poisson_background_target_mu"] = (
                args.poisson_background_target_mu
            )
            output.attrs["poisson_background_seed"] = args.poisson_background_seed
            output.attrs["poisson_background_preserve_zero_voxels"] = True
        background_original_ds = output.create_dataset(
            "background_mu_original", shape=(n,), dtype=np.float32
        )
        background_effective_ds = output.create_dataset(
            "background_mu_effective", shape=(n,), dtype=np.float32
        )
        background_added_ds = output.create_dataset(
            "background_mu_added", shape=(n,), dtype=np.float32
        )

        for event_index in range(n):
            volume = np.asarray(volumes[event_index])
            event_id = int(source["event_id"][event_index])
            original_background_mu = float(source["background_mu"][event_index])
            background_mu = original_background_mu
            added_background_mu = 0.0
            if args.poisson_background_target_mu is not None:
                event_seed = int(
                    np.random.SeedSequence(
                        [args.poisson_background_seed, event_id, 93563]
                    ).generate_state(1)[0]
                )
                volume, background_mu, added_background_mu = augment_volume_background(
                    volume,
                    original_background_mu,
                    args.poisson_background_target_mu,
                    seed=event_seed,
                )
            background_original_ds[event_index] = original_background_mu
            background_effective_ds[event_index] = background_mu
            background_added_ds[event_index] = added_background_mu
            score_grid = np.empty((len(y_positions), len(x_positions)), dtype=np.float32)
            slope_grid = np.full((*score_grid.shape, 2), np.nan, dtype=np.float32)
            deferred_regression = (
                regression_model is not model and args.regression_min_score is not None
            )
            with torch.no_grad():
                for locations, batch in crop_batches(
                    volume, background_mu, y_positions, x_positions,
                    args.crop_size, args.batch_size, normalization_mode,
                ):
                    batch = batch.to(device)
                    prediction = model(batch)
                    batch_scores = torch.sigmoid(prediction["presence_logit"]).cpu().numpy()
                    batch_slopes = None
                    if not deferred_regression:
                        regression_prediction = (
                            prediction if regression_model is model else regression_model(batch)
                        )
                        batch_slopes = regression_prediction["slope_xy"].cpu().numpy()
                    for batch_index, (row, col) in enumerate(locations):
                        score_grid[row, col] = batch_scores[batch_index]
                        if batch_slopes is not None:
                            slope_grid[row, col] = batch_slopes[batch_index]
                if deferred_regression:
                    selected_locations = list(zip(
                        *np.nonzero(score_grid >= args.regression_min_score), strict=True
                    ))
                    for first in range(0, len(selected_locations), args.batch_size):
                        locations = selected_locations[first : first + args.batch_size]
                        crops = np.stack([
                            normalize_scan_crop(
                                volume, background_mu, y_positions[row], x_positions[col],
                                args.crop_size, normalization_mode,
                            )
                            for row, col in locations
                        ])
                        batch = torch.from_numpy(crops[:, None]).float().to(device)
                        batch_slopes = regression_model(batch)["slope_xy"].cpu().numpy()
                        for batch_index, (row, col) in enumerate(locations):
                            slope_grid[row, col] = batch_slopes[batch_index]
            scores_ds[event_index] = score_grid
            slopes_ds[event_index] = slope_grid
            x_edges = np.asarray(source["x_edges_um"][event_index])
            y_edges = np.asarray(source["y_edges_um"][event_index])
            has_truth = "truth" in source
            if has_truth:
                truth_x_bin, truth_y_bin = truth_reference_bins(
                    source, event_index, x_edges, y_edges, plates,
                    args.truth_propagation_limit,
                )
            else:
                truth_x_bin = truth_y_bin = float("nan")
            components = connected_components(score_grid >= args.threshold)
            matched_candidates = 0
            for candidate_id, component in enumerate(components):
                representative = max(component, key=lambda location: float(score_grid[location]))
                row, col = representative
                slope = slope_grid[row, col]
                theta, phi = slopes_to_angles(slope[None])
                truth_matched = bool(
                    has_truth and cluster_contains_truth(
                        component, truth_x_bin, truth_y_bin,
                        x_positions, y_positions, args.crop_size,
                    )
                )
                matched_candidates += int(truth_matched)
                row_values = {
                    "event_index": event_index, "event_id": event_id,
                    "candidate_id": candidate_id, "grid_row": row, "grid_col": col,
                    "cluster_size": len(component),
                    "presence_score": float(score_grid[row, col]),
                    "slope_x": float(slope[0]), "slope_y": float(slope[1]),
                    "theta_pred_mrad": float(theta[0]), "phi_pred_rad": float(phi[0]),
                    "x_center_um": local_crop_center(x_edges, x_positions[col], args.crop_size),
                    "y_center_um": local_crop_center(y_edges, y_positions[row], args.crop_size),
                    "truth_matched": truth_matched,
                }
                for key, value in row_values.items():
                    candidates[key].append(value)
            event_report = {
                "event_id": event_id,
                "windows": int(score_grid.size),
                "windows_above_threshold": int((score_grid >= args.threshold).sum()),
                "candidate_clusters": len(components),
                "truth_matched_clusters": matched_candidates,
                "truth_detected": bool(matched_candidates > 0) if has_truth else None,
                "maximum_score": float(score_grid.max()),
            }
            for metadata_name in ("cell_id", "cell_x", "cell_y", "run_wall_brick"):
                if metadata_name in source:
                    event_report[metadata_name] = int(source[metadata_name][event_index])
            event_reports.append(event_report)

        group = output.create_group("candidates")
        integer_fields = {"event_index", "event_id", "candidate_id", "grid_row", "grid_col", "cluster_size"}
        for key, values in candidates.items():
            if key in integer_fields:
                dtype = np.int64
            elif key == "truth_matched":
                dtype = np.bool_
            else:
                dtype = np.float32
            group.create_dataset(key, data=np.asarray(values, dtype=dtype))

    report = {
        "input": str(args.input), "output": str(args.output),
        "crop_size": args.crop_size, "stride": args.stride,
        "grid_shape": [len(y_positions), len(x_positions)],
        "presence_threshold": args.threshold,
        "truth_propagation_limit": args.truth_propagation_limit,
        "classification_checkpoint": str(args.checkpoint),
        "regression_checkpoint": str(args.regression_checkpoint or args.checkpoint),
        "regression_min_score": args.regression_min_score,
        "device": str(device),
        "poisson_background_target_mu": args.poisson_background_target_mu,
        "poisson_background_seed": args.poisson_background_seed,
        "poisson_background_preserve_zero_voxels": True,
        "events": event_reports,
        "total_candidate_clusters": len(candidates["event_id"]),
    }
    args.output.with_suffix(".report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
