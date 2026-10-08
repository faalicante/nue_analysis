#!/usr/bin/env python3
"""Rebuild centred signal ROOT crops directly from 57x200x200 scan shards."""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import h5py
import numpy as np
import uproot

from scanning.root_xypseg_to_scan_hdf5 import crop_bounds


N_LAYERS = 57
PLATE_INDEX_OFFSET = 3
Z_BIN_SIZE_UM = 1350.0
REFERENCE_BRANCHES = (
    "background_mu", "zScale", "presence", "sample_type", "slope_x", "slope_y",
    "signal_event_id", "cell_id", "tag_cell_id",
)


def index_scan_events(paths: list[Path]) -> dict[int, tuple[Path, int]]:
    index: dict[int, tuple[Path, int]] = {}
    for path in paths:
        with h5py.File(path, "r") as source:
            for local_index, event_id in enumerate(np.asarray(source["event_id"], dtype=np.int64)):
                event = int(event_id)
                if event in index:
                    raise ValueError(f"duplicate event_id={event} in scan shards")
                index[event] = (path, local_index)
    return index


def projected_center_um(
    source: h5py.File,
    index: int,
    propagation_limit: int,
) -> tuple[float, float, int, int]:
    crop_plate = int(source["truth/plate"][index]) - PLATE_INDEX_OFFSET
    propagated = min(propagation_limit, N_LAYERS - crop_plate)
    if propagated < 0:
        raise ValueError(f"negative propagation for crop plate {crop_plate}")
    factor = Z_BIN_SIZE_UM * 0.5 * propagated / 1000.0
    center_x = (
        float(source["truth/xpos_um"][index])
        + float(source["truth/ltx_mrad_corrected"][index]) * factor
    )
    center_y = (
        float(source["truth/ypos_um"][index])
        + float(source["truth/lty_mrad_corrected"][index]) * factor
    )
    return center_x, center_y, crop_plate, propagated


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-glob", required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--propagation-limit", type=int, default=24)
    parser.add_argument("--crop-size", type=int, default=32)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--max-events", type=int)
    parser.add_argument(
        "--skip-missing-scan-events",
        action="store_true",
        help="Skip reference events absent from the scan shards instead of failing",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.crop_size <= 0 or args.crop_size % 2:
        raise ValueError("crop size must be a positive even integer")
    if args.propagation_limit <= 0:
        raise ValueError("propagation limit must be positive")
    if args.batch_size <= 0:
        raise ValueError("batch size must be positive")
    if args.output_root.exists() and not args.overwrite:
        raise FileExistsError(f"output exists: {args.output_root}")

    scan_paths = sorted(map(Path, glob.glob(args.scan_glob)))
    if not scan_paths:
        raise ValueError("scan glob matched no HDF5 files")
    scan_index = index_scan_events(scan_paths)

    with uproot.open(args.reference_root) as reference_file:
        tree = reference_file["samples"]
        event_ids = np.asarray(tree["signal_event_id"].array(library="np"), dtype=np.int64)
        metadata = {
            branch: np.asarray(tree[branch].array(library="np"))
            for branch in REFERENCE_BRANCHES
        }
    if args.max_events is not None:
        if args.max_events <= 0:
            raise ValueError("max-events must be positive")
        event_ids = event_ids[:args.max_events]
        metadata = {branch: values[:args.max_events] for branch, values in metadata.items()}
    missing = sorted(set(map(int, event_ids)) - set(scan_index))
    if missing and not args.skip_missing_scan_events:
        raise ValueError(f"{len(missing)} reference events absent from scan shards: {missing[:20]}")
    if missing:
        keep = np.asarray([int(event_id) in scan_index for event_id in event_ids], dtype=bool)
        event_ids = event_ids[keep]
        metadata = {branch: values[keep] for branch, values in metadata.items()}

    args.output_root.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output_root.with_suffix(".partial.root")
    if partial.exists():
        partial.unlink()
    branch_types: dict[str, object] = {
        "counts": np.dtype((np.int32, (N_LAYERS, args.crop_size, args.crop_size))),
    }
    branch_types.update({branch: values.dtype for branch, values in metadata.items()})
    handles: dict[Path, h5py.File] = {}
    diagnostics: list[dict[str, object]] = []
    try:
        with uproot.recreate(partial) as output_file:
            output_tree = output_file.mktree("samples", branch_types)
            for first in range(0, len(event_ids), args.batch_size):
                last = min(first + args.batch_size, len(event_ids))
                counts: list[np.ndarray] = []
                for output_index in range(first, last):
                    event_id = int(event_ids[output_index])
                    path, local_index = scan_index[event_id]
                    if path not in handles:
                        handles[path] = h5py.File(path, "r")
                    source = handles[path]
                    center_x, center_y, crop_plate, propagated = projected_center_um(
                        source, local_index, args.propagation_limit
                    )
                    x_start, x_stop, _ = crop_bounds(
                        np.asarray(source["x_edges_um"][local_index]), center_x, args.crop_size
                    )
                    y_start, y_stop, _ = crop_bounds(
                        np.asarray(source["y_edges_um"][local_index]), center_y, args.crop_size
                    )
                    crop = np.asarray(
                        source["volumes_raw"][local_index, :, y_start:y_stop, x_start:x_stop],
                        dtype=np.int32,
                    )
                    if crop.shape != (N_LAYERS, args.crop_size, args.crop_size):
                        raise ValueError(f"event {event_id}: unexpected crop shape {crop.shape}")
                    reference_mu = float(metadata["background_mu"][output_index])
                    scan_mu = float(source["background_mu"][local_index])
                    if not np.isclose(reference_mu, scan_mu, rtol=0.0, atol=1e-5):
                        raise ValueError(
                            f"event {event_id}: background_mu differs ({reference_mu} vs {scan_mu})"
                        )
                    counts.append(crop)
                    if len(diagnostics) < 20:
                        diagnostics.append({
                            "event_id": event_id,
                            "source_shard": str(path),
                            "source_index": local_index,
                            "crop_plate": crop_plate,
                            "propagated_plates": propagated,
                            "center_xy_um": [center_x, center_y],
                            "start_xy_bin_zero_based": [x_start, y_start],
                            "crop_sum": int(crop.sum()),
                        })
                batch = {"counts": np.stack(counts)}
                batch.update({branch: values[first:last] for branch, values in metadata.items()})
                output_tree.extend(batch)
                print(json.dumps({"written": last, "total": len(event_ids)}), flush=True)
    finally:
        for handle in handles.values():
            handle.close()

    partial.replace(args.output_root)
    report = {
        "output_root": str(args.output_root),
        "reference_root": str(args.reference_root),
        "scan_glob": args.scan_glob,
        "scan_shards": len(scan_paths),
        "events": len(event_ids),
        "skipped_reference_events_absent_from_scan": missing,
        "crop_shape": [N_LAYERS, args.crop_size, args.crop_size],
        "crop_center_formula": (
            "p0=truth_plate-3; n=min(propagation_limit,57-p0); "
            "center=vertex+slope_corrected*1350um*0.5*n/1000"
        ),
        "propagation_limit": args.propagation_limit,
        "metadata_source": "reference ROOT; only counts are rebuilt from scan shards",
        "diagnostics_first_events": diagnostics,
        "size_bytes": args.output_root.stat().st_size,
    }
    report_path = args.output_root.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
