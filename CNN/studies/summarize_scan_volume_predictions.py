#!/usr/bin/env python3
"""Summarize event-level scan efficiency and candidate multiplicity versus threshold."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from training.cnn21d.metrics import slopes_to_angles
from scanning.scan_cnn21d_volumes import cluster_contains_truth, connected_components


ANGLE_BINS = (
    ("5_to_10", 5.0, 10.0),
    ("10_to_20", 10.0, 20.0),
    ("20_to_50", 20.0, 50.0),
    ("50_to_100", 50.0, 100.0),
    ("ge_100", 100.0, np.inf),
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--volumes", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--thresholds", type=float, nargs="+",
        default=[0.5, 0.9, 0.95, 0.98, 0.99, 0.995],
    )
    parser.add_argument(
        "--truth-propagation-limit", type=int,
        help=(
            "recompute the truth-matching point with "
            "n=min(limit,57-(truth_plate-3)); otherwise use the point stored in the HDF5"
        ),
    )
    args = parser.parse_args()

    with h5py.File(args.volumes, "r") as source, h5py.File(args.predictions, "r") as prediction:
        event_ids = source["event_id"][:]
        truth_theta = 1000.0 * np.arctan(
            np.hypot(
                source["truth/ltx_mrad_corrected"][:],
                source["truth/lty_mrad_corrected"][:],
            ) / 1000.0
        )
        scores = prediction["window_presence_score"][:]
        slopes = prediction["window_slope_xy"][:]
        x_positions = prediction["grid_x_start_bin"][:].tolist()
        y_positions = prediction["grid_y_start_bin"][:].tolist()
        crop_size = int(prediction.attrs["crop_size"])
        if args.truth_propagation_limit is None:
            truth_x = source["truth/propagated_center_x_bin"][:]
            truth_y = source["truth/propagated_center_y_bin"][:]
            truth_matching = "stored"
        else:
            if args.truth_propagation_limit <= 0:
                raise ValueError("--truth-propagation-limit must be positive")
            crop_plate = source["truth/plate"][:].astype(np.int64) - 3
            propagated = np.minimum(args.truth_propagation_limit, 57 - crop_plate)
            factor = 1350.0 * 0.5 * propagated / 1000.0
            center_x_um = (
                source["truth/xpos_um"][:]
                + source["truth/ltx_mrad_corrected"][:] * factor
            )
            center_y_um = (
                source["truth/ypos_um"][:]
                + source["truth/lty_mrad_corrected"][:] * factor
            )
            truth_x = np.asarray([
                np.interp(value, 0.5 * (edges[:-1] + edges[1:]), np.arange(len(edges) - 1))
                for value, edges in zip(center_x_um, source["x_edges_um"][:], strict=True)
            ])
            truth_y = np.asarray([
                np.interp(value, 0.5 * (edges[:-1] + edges[1:]), np.arange(len(edges) - 1))
                for value, edges in zip(center_y_um, source["y_edges_um"][:], strict=True)
            ])
            truth_matching = f"recomputed_propagation_limit_{args.truth_propagation_limit}"

        summaries: list[dict[str, object]] = []
        for threshold in args.thresholds:
            detected = np.zeros(len(event_ids), dtype=bool)
            theta_prediction = np.full(len(event_ids), np.nan, dtype=np.float64)
            total_clusters = 0
            unmatched_clusters = 0
            events_with_unmatched = 0
            for index in range(len(event_ids)):
                components = connected_components(scores[index] >= threshold)
                total_clusters += len(components)
                event_has_unmatched = False
                best_matched_score = -np.inf
                for component in components:
                    representative = max(component, key=lambda rc: float(scores[index][rc]))
                    matched = cluster_contains_truth(
                        component, float(truth_x[index]), float(truth_y[index]),
                        x_positions, y_positions, crop_size,
                    )
                    if matched:
                        detected[index] = True
                        representative_score = float(scores[index][representative])
                        if representative_score > best_matched_score:
                            best_matched_score = representative_score
                            slope = slopes[index][representative]
                            theta_prediction[index] = float(slopes_to_angles(slope[None])[0][0])
                    else:
                        unmatched_clusters += 1
                        event_has_unmatched = True
                events_with_unmatched += int(event_has_unmatched)

            by_angle: dict[str, object] = {}
            for name, low, high in ANGLE_BINS:
                mask = (truth_theta >= low) & (truth_theta < high)
                regression_mask = mask & detected
                residual = theta_prediction[regression_mask] - truth_theta[regression_mask]
                by_angle[name] = {
                    "events": int(mask.sum()),
                    "detected": int((mask & detected).sum()),
                    "efficiency": float(detected[mask].mean()) if mask.any() else None,
                    "theta_mae_mrad_detected": (
                        float(np.abs(residual).mean()) if residual.size else None
                    ),
                }
            summaries.append({
                "threshold": threshold,
                "events": len(event_ids),
                "detected_events": int(detected.sum()),
                "total_candidate_clusters": total_clusters,
                "unmatched_candidate_clusters": unmatched_clusters,
                "events_with_unmatched_candidate": events_with_unmatched,
                "angle_bins": by_angle,
            })

    result = {
        "volumes": str(args.volumes),
        "predictions": str(args.predictions),
        "truth_matching": truth_matching,
        "operating_points": summaries,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
