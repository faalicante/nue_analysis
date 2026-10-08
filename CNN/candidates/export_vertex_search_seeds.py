#!/usr/bin/env python3
"""Export inclined-shower candidates as backward tracks for vertex matching.

The input can be a data/background GIF manifest or a signal XZ/YZ manifest.
For the latter, pass ``--volumes`` so event ids can be mapped back to their
source HDF5 entries. The MC residual table calibrates the longitudinal and
transverse search windows.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import h5py
import numpy as np


Z_STEP_UM = 1350.0
NOMINAL_SHOWER_DELAY_PLATES = 3
FIRST_XYPSEG_LAYER = 1
LAST_XYPSEG_LAYER = 57


def decode(value: object) -> object:
    return value.decode("utf-8") if isinstance(value, bytes) else value


def quantile(values: np.ndarray, probability: float) -> float:
    return float(np.quantile(values, probability, method="linear"))


def vertex_interval(
    visible_layer: int,
    residual_low: float,
    residual_high: float,
) -> tuple[int, int, int, int]:
    """Return raw and detector-clipped vertex-layer bounds.

    residual = reconstructed visible layer - nominal visible layer, while the
    nominal visible layer is vertex + 3.  A physical vertex is constrained to
    be at least three plates upstream of the reconstructed visible component.
    """
    raw_low = math.floor(
        visible_layer - NOMINAL_SHOWER_DELAY_PLATES - residual_high
    )
    raw_high = math.ceil(
        visible_layer - NOMINAL_SHOWER_DELAY_PLATES - residual_low
    )
    raw_high = min(raw_high, visible_layer - NOMINAL_SHOWER_DELAY_PLATES)
    low = max(FIRST_XYPSEG_LAYER, raw_low)
    high = min(LAST_XYPSEG_LAYER, raw_high)
    return raw_low, raw_high, low, high


def project_xy(
    start_x_um: float,
    start_y_um: float,
    tx_mrad: float,
    ty_mrad: float,
    start_layer: int,
    target_layer: int,
) -> tuple[float, float]:
    delta_z_um = (target_layer - start_layer) * Z_STEP_UM
    return (
        start_x_um + tx_mrad * delta_z_um / 1000.0,
        start_y_um + ty_mrad * delta_z_um / 1000.0,
    )


def source_metadata(volume_path: str, event_index: int) -> dict[str, object]:
    with h5py.File(volume_path, "r") as source:
        x_edges = np.asarray(source["x_edges_um"][event_index], dtype=float)
        y_edges = np.asarray(source["y_edges_um"][event_index], dtype=float)
        return {
            "source_root_file": str(decode(source["source_file"][event_index])),
            "xy_bin_width_x_um": float(np.median(np.diff(x_edges))),
            "xy_bin_width_y_um": float(np.median(np.diff(y_edges))),
        }


def event_lookup_from_volume(volume_path: Path) -> dict[int, tuple[str, int]]:
    with h5py.File(volume_path, "r") as source:
        event_ids = np.asarray(source["event_id"], dtype=np.int64)
    lookup: dict[int, tuple[str, int]] = {}
    for index, event_id in enumerate(event_ids):
        event_id = int(event_id)
        if event_id in lookup:
            raise ValueError(f"event_id {event_id} appears more than once in {volume_path}")
        lookup[event_id] = (str(volume_path), index)
    return lookup


def candidate_source_location(
    candidate: dict[str, object], event_lookup: dict[int, tuple[str, int]],
) -> tuple[str, int]:
    """Resolve a candidate's source from its own fields or an event-id lookup."""
    if "volume_path" in candidate and "event_index" in candidate:
        return str(candidate["volume_path"]), int(candidate["event_index"])
    event_id = int(candidate["event_id"])
    if event_id not in event_lookup:
        raise ValueError(
            "manifest has no volume_path/event_index; pass --volumes containing "
            f"event_id {event_id}"
        )
    return event_lookup[event_id]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--volumes", type=Path,
        help="Signal scan HDF5 used to resolve source entries missing from a signal manifest.",
    )
    parser.add_argument("--mc-residuals", type=Path, required=True)
    parser.add_argument("--mc-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-slope-disagreement-mrad", type=float, default=5.0)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    candidates = manifest["candidates"]
    event_lookup = event_lookup_from_volume(args.volumes) if args.volumes else {}
    with args.mc_residuals.open(newline="", encoding="utf-8") as handle:
        residual_rows = list(csv.DictReader(handle))
    residuals = np.asarray(
        [float(row["layer_residual"]) for row in residual_rows], dtype=float
    )
    mc_summary = json.loads(args.mc_summary.read_text(encoding="utf-8"))
    position = mc_summary["position_error_at_nominal_start_um"]

    q05, q16, q50, q84, q95 = (
        quantile(residuals, q) for q in (0.05, 0.16, 0.50, 0.84, 0.95)
    )
    radius_core_um = float(position["q68"])
    radius_broad_um = float(position["q90"])

    seed_rows: list[dict[str, object]] = []
    plate_rows: list[dict[str, object]] = []
    metadata_cache: dict[tuple[str, int], dict[str, object]] = {}
    for candidate in candidates:
        volume_path, event_index = candidate_source_location(candidate, event_lookup)
        key = (volume_path, event_index)
        if key not in metadata_cache:
            metadata_cache[key] = source_metadata(volume_path, event_index)
        source = metadata_cache[key]

        start_layer = int(candidate["start_layer"])
        start_x = float(candidate["start_x_um"])
        start_y = float(candidate["start_y_um"])
        cnn_tx = float(
            candidate["tx_pred_mrad"] if "tx_pred_mrad" in candidate else candidate["tx_cnn_mrad"]
        )
        cnn_ty = float(
            candidate["ty_pred_mrad"] if "ty_pred_mrad" in candidate else candidate["ty_cnn_mrad"]
        )
        if (
            candidate.get("selected_component_slope_x") is not None
            and candidate.get("selected_component_slope_y") is not None
        ):
            component_tx = (
                1000.0 * float(candidate["selected_component_slope_x"])
                * float(source["xy_bin_width_x_um"]) / Z_STEP_UM
            )
            component_ty = (
                1000.0 * float(candidate["selected_component_slope_y"])
                * float(source["xy_bin_width_y_um"]) / Z_STEP_UM
            )
        else:
            # XZ/YZ manifests already store the fitted slope in physical mrad.
            component_tx = float(candidate.get("tx_fit_mrad", float("nan")))
            component_ty = float(candidate.get("ty_fit_mrad", float("nan")))
        slope_disagreement = float(math.hypot(component_tx - cnn_tx, component_ty - cnn_ty))
        use_component = (
            candidate.get("start_quality") == "coherent"
            and np.isfinite(component_tx)
            and np.isfinite(component_ty)
            and slope_disagreement <= args.max_slope_disagreement_mrad
        )
        recommended_tx = component_tx if use_component else cnn_tx
        recommended_ty = component_ty if use_component else cnn_ty
        recommended_source = "fitted_hot_component" if use_component else "cnn"

        core_raw_lo, core_raw_hi, core_lo, core_hi = vertex_interval(
            start_layer, q16, q84
        )
        broad_raw_lo, broad_raw_hi, broad_lo, broad_hi = vertex_interval(
            start_layer, q05, q95
        )
        central_raw = int(round(start_layer - NOMINAL_SHOWER_DELAY_PLATES - q50))
        central = min(LAST_XYPSEG_LAYER, max(FIRST_XYPSEG_LAYER, central_raw))
        central_x, central_y = project_xy(
            start_x, start_y, recommended_tx, recommended_ty, start_layer, central
        )

        seed = {
            "rank": int(candidate["rank"]),
            "brick_name": str(candidate.get("brick_name", "signal_mc_validation")),
            "run_wall_brick": int(candidate.get("run_wall_brick", -1)),
            "event_id": int(candidate["event_id"]),
            "cell_id": int(candidate["cell_id"]),
            "cell_x": int(candidate["cell_x"]),
            "cell_y": int(candidate["cell_y"]),
            "candidate_id": int(candidate["candidate_id"]),
            "source_root_file": source["source_root_file"],
            "presence_score": float(candidate["presence_score"]),
            "theta_cnn_mrad": float(
                candidate["theta_pred_mrad"]
                if "theta_pred_mrad" in candidate else candidate["theta_cnn_mrad"]
            ),
            "cnn_tx_mrad": cnn_tx,
            "cnn_ty_mrad": cnn_ty,
            "component_tx_mrad": component_tx,
            "component_ty_mrad": component_ty,
            "component_theta_mrad": float(
                1000.0 * math.atan(math.hypot(component_tx, component_ty) / 1000.0)
            ),
            "cnn_component_slope_disagreement_mrad": slope_disagreement,
            "component_slope_compatible": use_component,
            "recommended_slope_source": recommended_source,
            "recommended_tx_mrad": recommended_tx,
            "recommended_ty_mrad": recommended_ty,
            "visible_start_layer": start_layer,
            "visible_start_root_z_bin": int(candidate.get("start_root_z_bin", start_layer + 1)),
            "visible_start_x_um": start_x,
            "visible_start_y_um": start_y,
            "vertex_central_layer": central,
            "vertex_central_root_z_bin": central + 1,
            "vertex_central_x_um": central_x,
            "vertex_central_y_um": central_y,
            "vertex_core_layer_first": core_lo if core_lo <= core_hi else "",
            "vertex_core_layer_last": core_hi if core_lo <= core_hi else "",
            "vertex_broad_layer_first": broad_lo if broad_lo <= broad_hi else "",
            "vertex_broad_layer_last": broad_hi if broad_lo <= broad_hi else "",
            "vertex_core_radius_um": radius_core_um,
            "vertex_broad_radius_um": radius_broad_um,
            "upstream_search_truncated": bool(broad_raw_lo < FIRST_XYPSEG_LAYER),
            "core_raw_layer_first": core_raw_lo,
            "core_raw_layer_last": core_raw_hi,
            "broad_raw_layer_first": broad_raw_lo,
            "broad_raw_layer_last": broad_raw_hi,
            "start_quality": candidate["start_quality"],
            "start_component_voxels": candidate.get("start_component_voxels", ""),
            "start_component_layers": candidate.get("start_component_layers", ""),
            "start_component_charge": candidate.get("start_component_charge", ""),
            "gif": candidate.get("gif", candidate.get("image", "")),
        }
        seed_rows.append(seed)

        if broad_lo <= broad_hi:
            for layer in range(broad_lo, broad_hi + 1):
                rec_x, rec_y = project_xy(
                    start_x, start_y, recommended_tx, recommended_ty, start_layer, layer
                )
                cnn_x, cnn_y = project_xy(
                    start_x, start_y, cnn_tx, cnn_ty, start_layer, layer
                )
                component_x, component_y = project_xy(
                    start_x, start_y, component_tx, component_ty, start_layer, layer
                )
                plate_rows.append({
                    "rank": seed["rank"],
                    "brick_name": seed["brick_name"],
                    "event_id": seed["event_id"],
                    "cell_id": seed["cell_id"],
                    "cell_x": seed["cell_x"],
                    "cell_y": seed["cell_y"],
                    "candidate_id": seed["candidate_id"],
                    "source_root_file": seed["source_root_file"],
                    "vertex_layer": layer,
                    "vertex_root_z_bin": layer + 1,
                    "in_core_interval": bool(core_lo <= layer <= core_hi),
                    "recommended_x_um": rec_x,
                    "recommended_y_um": rec_y,
                    "recommended_slope_source": recommended_source,
                    "cnn_x_um": cnn_x,
                    "cnn_y_um": cnn_y,
                    "component_x_um": component_x,
                    "component_y_um": component_y,
                    "core_radius_um": radius_core_um,
                    "broad_radius_um": radius_broad_um,
                })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    seeds_path = args.output_dir / "vertex_search_seeds.csv"
    plates_path = args.output_dir / "vertex_search_by_plate.csv"
    with seeds_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(seed_rows[0]))
        writer.writeheader()
        writer.writerows(seed_rows)
    with plates_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(plate_rows[0]))
        writer.writeheader()
        writer.writerows(plate_rows)

    calibration = {
        "candidate_count": len(seed_rows),
        "plate_search_rows": len(plate_rows),
        "layer_convention": {
            "xypseg_layer": "1..57",
            "root_z_bin": "xypseg_layer + 1 (2..58)",
            "one_xypseg_layer_per_physical_plate": True,
        },
        "vertex_model": (
            "visible_start_layer = vertex_layer + 3 + MC_residual; all output "
            "intervals enforce vertex_layer <= visible_start_layer - 3"
        ),
        "mc_layer_residual_quantiles": {
            "q05": q05, "q16": q16, "q50": q50, "q84": q84, "q95": q95
        },
        "spatial_radius_um": {
            "core_q68": radius_core_um,
            "broad_q90": radius_broad_um,
        },
        "recommended_slope": (
            "fitted Poisson-hot component only when its vector differs from the "
            f"CNN by <= {args.max_slope_disagreement_mrad:g} mrad; otherwise CNN"
        ),
        "outputs": {
            "candidate_seeds": str(seeds_path.resolve()),
            "one_row_per_candidate_plate": str(plates_path.resolve()),
        },
    }
    calibration_path = args.output_dir / "vertex_search_calibration.json"
    calibration_path.write_text(json.dumps(calibration, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(calibration, indent=2))


if __name__ == "__main__":
    main()
