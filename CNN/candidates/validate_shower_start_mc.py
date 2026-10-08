#!/usr/bin/env python3
"""Validate the diagnostic shower-start estimator against signal MC metadata."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np

from visualization.make_background_candidate_gifs import centered_window_start, estimate_shower_start


Z_STEP_UM = 1350.0


def metric_summary(values: np.ndarray) -> dict[str, float]:
    absolute = np.abs(values)
    return {
        "bias": float(np.mean(values)),
        "median": float(np.median(values)),
        "mae": float(np.mean(absolute)),
        "absolute_q68": float(np.quantile(absolute, 0.68)),
        "absolute_q90": float(np.quantile(absolute, 0.90)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--volumes", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--theta-min-mrad", type=float, default=10.0)
    parser.add_argument("--display-crop-size", type=int, default=40)
    parser.add_argument("--shower-delay-plates", type=float, default=3.0)
    parser.add_argument("--xypseg-layers-per-plate", type=int, default=1)
    parser.add_argument("--max-slope-disagreement-mrad", type=float, default=5.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    with h5py.File(args.volumes, "r") as source, h5py.File(args.predictions, "r") as pred:
        group = pred["candidates"]
        x_positions = np.asarray(pred["grid_x_start_bin"], dtype=int)
        y_positions = np.asarray(pred["grid_y_start_bin"], dtype=int)
        cnn_crop_size = int(pred.attrs["crop_size"])
        best: dict[int, int] = {}
        for candidate_index in range(len(group["event_id"])):
            if not bool(group["truth_matched"][candidate_index]):
                continue
            event_index = int(group["event_index"][candidate_index])
            previous = best.get(event_index)
            if previous is None or group["presence_score"][candidate_index] > group["presence_score"][previous]:
                best[event_index] = candidate_index

        delay_layers = int(round(args.shower_delay_plates * args.xypseg_layers_per_plate))
        for event_index, candidate_index in best.items():
            tx = float(source["truth/ltx_mrad_corrected"][event_index])
            ty = float(source["truth/lty_mrad_corrected"][event_index])
            theta = float(1000.0 * np.arctan(np.hypot(tx, ty) / 1000.0))
            if theta <= args.theta_min_mrad:
                continue
            x0 = int(x_positions[int(group["grid_col"][candidate_index])])
            y0 = int(y_positions[int(group["grid_row"][candidate_index])])
            full = source["volumes_raw"].shape[-1]
            display_x0 = centered_window_start(
                x0, cnn_crop_size, args.display_crop_size, full
            )
            display_y0 = centered_window_start(
                y0, cnn_crop_size, args.display_crop_size, full
            )
            raw = np.asarray(source["volumes_raw"][
                event_index, :,
                display_y0:display_y0 + args.display_crop_size,
                display_x0:display_x0 + args.display_crop_size,
            ])
            x_edges = np.asarray(source["x_edges_um"][event_index])
            y_edges = np.asarray(source["y_edges_um"][event_index])
            x_centers = 0.5 * (
                x_edges[display_x0:display_x0 + args.display_crop_size]
                + x_edges[display_x0 + 1:display_x0 + args.display_crop_size + 1]
            )
            y_centers = 0.5 * (
                y_edges[display_y0:display_y0 + args.display_crop_size]
                + y_edges[display_y0 + 1:display_y0 + args.display_crop_size + 1]
            )
            estimate = estimate_shower_start(
                raw,
                float(source["background_mu"][event_index]),
                x_centers,
                y_centers,
                predicted_slope_xy=np.asarray([
                    group["slope_x"][candidate_index], group["slope_y"][candidate_index]
                ]),
                candidate_center_xy_um=np.asarray([
                    group["x_center_um"][candidate_index],
                    group["y_center_um"][candidate_index],
                ]),
            )
            # MC stores the vertex plate. In XYPseg indexing the vertex is
            # truth_plate-3 and one XYPseg layer corresponds to one plate, so a
            # shower delayed by 3 plates starts at truth_plate-3+3=truth_plate.
            expected_layer = int(source["truth/plate"][event_index]) - 3 + delay_layers
            vertex_layer = int(source["truth/plate"][event_index]) - 3
            x_bin_width_um = float(np.median(np.diff(x_edges)))
            y_bin_width_um = float(np.median(np.diff(y_edges)))
            cnn_tx = (
                1000.0 * float(group["slope_x"][candidate_index])
                * x_bin_width_um / Z_STEP_UM
            )
            cnn_ty = (
                1000.0 * float(group["slope_y"][candidate_index])
                * y_bin_width_um / Z_STEP_UM
            )
            component_tx = (
                1000.0 * float(estimate["selected_component_slope_x"])
                * x_bin_width_um / Z_STEP_UM
            )
            component_ty = (
                1000.0 * float(estimate["selected_component_slope_y"])
                * y_bin_width_um / Z_STEP_UM
            )
            slope_disagreement = float(np.hypot(
                component_tx - cnn_tx, component_ty - cnn_ty
            ))
            use_component = (
                estimate["start_quality"] == "coherent"
                and np.isfinite(component_tx)
                and np.isfinite(component_ty)
                and slope_disagreement <= args.max_slope_disagreement_mrad
            )
            recommended_tx = component_tx if use_component else cnn_tx
            recommended_ty = component_ty if use_component else cnn_ty
            recommended_source = "fitted_hot_component" if use_component else "cnn"
            projected_vertex_x = float(estimate["start_x_um"]) + (
                recommended_tx * Z_STEP_UM
                * (vertex_layer - int(estimate["start_layer"])) / 1000.0
            )
            projected_vertex_y = float(estimate["start_y_um"]) + (
                recommended_ty * Z_STEP_UM
                * (vertex_layer - int(estimate["start_layer"])) / 1000.0
            )
            vertex_projection_dx = projected_vertex_x - float(source["truth/xpos_um"][event_index])
            vertex_projection_dy = projected_vertex_y - float(source["truth/ypos_um"][event_index])
            distance_um = args.shower_delay_plates * Z_STEP_UM
            expected_x = float(source["truth/xpos_um"][event_index]) + tx * distance_um / 1000.0
            expected_y = float(source["truth/ypos_um"][event_index]) + ty * distance_um / 1000.0
            dx = float(estimate["start_x_um"]) - expected_x
            dy = float(estimate["start_y_um"]) - expected_y
            reconstructed_delta_layers = int(estimate["start_layer"]) - vertex_layer
            track_x_at_reconstructed_layer = (
                float(source["truth/xpos_um"][event_index])
                + tx * Z_STEP_UM * reconstructed_delta_layers / 1000.0
            )
            track_y_at_reconstructed_layer = (
                float(source["truth/ypos_um"][event_index])
                + ty * Z_STEP_UM * reconstructed_delta_layers / 1000.0
            )
            transverse_track_error = float(np.hypot(
                float(estimate["start_x_um"]) - track_x_at_reconstructed_layer,
                float(estimate["start_y_um"]) - track_y_at_reconstructed_layer,
            ))
            rows.append({
                "event_id": int(source["event_id"][event_index]),
                "presence_score": float(group["presence_score"][candidate_index]),
                "theta_true_mrad": theta,
                "expected_start_layer": expected_layer,
                "reconstructed_start_layer": int(estimate["start_layer"]),
                "layer_residual": int(estimate["start_layer"]) - expected_layer,
                "expected_start_x_um": expected_x,
                "expected_start_y_um": expected_y,
                "reconstructed_start_x_um": float(estimate["start_x_um"]),
                "reconstructed_start_y_um": float(estimate["start_y_um"]),
                "position_residual_x_um": dx,
                "position_residual_y_um": dy,
                "position_error_um": float(np.hypot(dx, dy)),
                "truth_track_x_at_reconstructed_layer_um": track_x_at_reconstructed_layer,
                "truth_track_y_at_reconstructed_layer_um": track_y_at_reconstructed_layer,
                "transverse_track_error_um": transverse_track_error,
                "cnn_tx_mrad": cnn_tx,
                "cnn_ty_mrad": cnn_ty,
                "component_tx_mrad": component_tx,
                "component_ty_mrad": component_ty,
                "cnn_component_slope_disagreement_mrad": slope_disagreement,
                "recommended_slope_source": recommended_source,
                "recommended_tx_mrad": recommended_tx,
                "recommended_ty_mrad": recommended_ty,
                "projected_vertex_x_um": projected_vertex_x,
                "projected_vertex_y_um": projected_vertex_y,
                "vertex_projection_residual_x_um": vertex_projection_dx,
                "vertex_projection_residual_y_um": vertex_projection_dy,
                "vertex_projection_error_um": float(np.hypot(
                    vertex_projection_dx, vertex_projection_dy
                )),
                **estimate,
            })

    layer_residual = np.asarray([row["layer_residual"] for row in rows], dtype=float)
    position_error = np.asarray([row["position_error_um"] for row in rows], dtype=float)
    coherent = np.asarray([row["start_quality"] == "coherent" for row in rows])
    transverse_error = np.asarray(
        [row["transverse_track_error_um"] for row in rows], dtype=float
    )
    vertex_projection_error = np.asarray(
        [row["vertex_projection_error_um"] for row in rows], dtype=float
    )
    layer_within_3 = np.abs(layer_residual) <= 3

    def position_summary(values: np.ndarray) -> dict[str, float]:
        return {
            "median": float(np.median(values)),
            "mean": float(np.mean(values)),
            "q68": float(np.quantile(values, 0.68)),
            "q90": float(np.quantile(values, 0.90)),
            "within_100_um": float(np.mean(values <= 100.0)),
            "within_250_um": float(np.mean(values <= 250.0)),
            "within_500_um": float(np.mean(values <= 500.0)),
        }
    summary = {
        "volumes": str(args.volumes),
        "predictions": str(args.predictions),
        "selection": f"best truth-matched candidate and theta_true > {args.theta_min_mrad} mrad",
        "reference_definition": {
            "vertex_layer": "truth_plate - 3",
            "shower_delay_physical_plates": args.shower_delay_plates,
            "xypseg_layers_per_plate": args.xypseg_layers_per_plate,
            "expected_shower_start_layer": "truth_plate - 3 + delay_layers",
        },
        "events": len(rows),
        "coherent_events": int(coherent.sum()),
        "layer_residual": metric_summary(layer_residual),
        "layer_accuracy": {
            "exact": float(np.mean(layer_residual == 0)),
            "within_1": float(np.mean(np.abs(layer_residual) <= 1)),
            "within_2": float(np.mean(np.abs(layer_residual) <= 2)),
            "within_3": float(np.mean(np.abs(layer_residual) <= 3)),
        },
        "position_error_at_nominal_start_um": position_summary(position_error),
        "transverse_error_to_truth_track_at_reconstructed_layer_um": position_summary(
            transverse_error
        ),
        "transverse_error_when_layer_within_3_um": {
            "events": int(layer_within_3.sum()),
            **position_summary(transverse_error[layer_within_3]),
        },
        "vertex_projection_to_original_mc_vertex_um": {
            **position_summary(vertex_projection_error),
        },
        "recommended_slope_source_counts": {
            source: int(sum(row["recommended_slope_source"] == source for row in rows))
            for source in ("cnn", "fitted_hot_component")
        },
    }
    csv_path = args.output_dir / "mc_shower_start_residuals.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
