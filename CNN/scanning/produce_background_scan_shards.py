#!/usr/bin/env python3
"""Produce restartable HDF5 shards for full-cell background CNN scans."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scanning.root_background_xypseg_to_scan_hdf5 import (
    BackgroundCell,
    DEFAULT_CELL_XMIN_UM,
    DEFAULT_CELL_YMIN_UM,
    render_background_root_path,
    write_background_hdf5,
    write_report,
)


DEFAULT_ROOT_TEMPLATE = (
    "/eos/experiment/sndlhc/users/dancc/FEDRA/"
    "muon_Euniform_RUN1_FLUKA25/b000021/cell_reco/"
    "cell_{x}0_{y}0/b000021/b000021.0.{x}.{y}.trk.root"
)
KNOWN_BAD_CELL_IDS = {58}
SCRIPT_VERSION = "2026-09-09-2"


def chunked(values: list[BackgroundCell], size: int) -> list[list[BackgroundCell]]:
    if size <= 0:
        raise ValueError("shard size must be positive")
    return [values[first:first + size] for first in range(0, len(values), size)]


def selected_cells(
    x_min: int,
    x_max: int,
    y_min: int,
    y_max: int,
    include_known_bad: bool,
    *,
    cell_xmin_um: float = DEFAULT_CELL_XMIN_UM,
    cell_ymin_um: float = DEFAULT_CELL_YMIN_UM,
) -> list[BackgroundCell]:
    if not (1 <= x_min <= x_max <= 18 and 1 <= y_min <= y_max <= 18):
        raise ValueError("x/y ranges must satisfy 1 <= min <= max <= 18")
    cells = [
        BackgroundCell(
            x=x,
            y=y,
            grid_xmin_um=cell_xmin_um,
            grid_ymin_um=cell_ymin_um,
        )
        for y in range(y_min, y_max + 1)
        for x in range(x_min, x_max + 1)
    ]
    if not include_known_bad:
        cells = [cell for cell in cells if cell.cell_id not in KNOWN_BAD_CELL_IDS]
    return cells


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-template", default=DEFAULT_ROOT_TEMPLATE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-size", type=int, default=25)
    parser.add_argument("--x-min", type=int, default=1)
    parser.add_argument("--x-max", type=int, default=18)
    parser.add_argument("--y-min", type=int, default=1)
    parser.add_argument("--y-max", type=int, default=18)
    parser.add_argument(
        "--cell-xmin-um",
        type=float,
        default=DEFAULT_CELL_XMIN_UM,
        help="Nominal lower x edge of cell (1,1), in micrometres",
    )
    parser.add_argument(
        "--cell-ymin-um",
        type=float,
        default=DEFAULT_CELL_YMIN_UM,
        help="Nominal lower y edge of cell (1,1), in micrometres",
    )
    parser.add_argument("--max-cells", type=int)
    parser.add_argument(
        "--include-known-bad-cell",
        action="store_true",
        help="Also process cell_id=58 (x=5,y=4), skipped by crop_bw.C for missing XYPseg data",
    )
    parser.add_argument(
        "--skip-errors",
        action="store_true",
        help="Record shards that fail during extraction and continue; default is fail-fast",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cells = selected_cells(
        args.x_min,
        args.x_max,
        args.y_min,
        args.y_max,
        args.include_known_bad_cell,
        cell_xmin_um=args.cell_xmin_um,
        cell_ymin_um=args.cell_ymin_um,
    )
    if args.max_cells is not None:
        if args.max_cells <= 0:
            raise ValueError("max-cells must be positive")
        cells = cells[:args.max_cells]
    if not cells:
        raise ValueError("cell selection is empty")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "manifest.json"
    manifest: dict[str, object] = {
        "producer_version": SCRIPT_VERSION,
        "root_template": args.root_template,
        "cell_xmin_um": args.cell_xmin_um,
        "cell_ymin_um": args.cell_ymin_um,
        "cell_ranges": {
            "x": [args.x_min, args.x_max],
            "y": [args.y_min, args.y_max],
        },
        "known_bad_cell_ids_excluded": (
            [] if args.include_known_bad_cell else sorted(KNOWN_BAD_CELL_IDS)
        ),
        "selected_cells": len(cells),
        "shard_size": args.shard_size,
        "shards": [],
    }

    for shard_index, shard_cells in enumerate(chunked(cells, args.shard_size)):
        first = shard_cells[0]
        last = shard_cells[-1]
        output = args.output_dir / (
            f"background_scan_{shard_index:04d}_"
            f"cell{first.cell_id:03d}-{last.cell_id:03d}.h5"
        )
        row: dict[str, object] = {
            "index": shard_index,
            "output": str(output),
            "cell_id_first": first.cell_id,
            "cell_id_last": last.cell_id,
            "requested_cell_count": len(shard_cells),
        }
        if output.exists() and not args.overwrite:
            row["status"] = "skipped_existing"
            row["cell_count"] = len(shard_cells)
        elif args.dry_run:
            row["status"] = "dry_run"
            row["cell_count"] = len(shard_cells)
            row["first_root"] = render_background_root_path(args.root_template, first)
        else:
            partial = output.with_suffix(".partial.h5")
            if partial.exists():
                partial.unlink()
            try:
                cells_and_paths = [
                    (cell, render_background_root_path(args.root_template, cell))
                    for cell in shard_cells
                ]
                report = write_background_hdf5(cells_and_paths, partial)
            except Exception as error:
                if partial.exists():
                    partial.unlink()
                failure = {
                    "scope": "shard",
                    "cell_id_first": first.cell_id,
                    "cell_id_last": last.cell_id,
                    "error_type": type(error).__name__,
                    "reason": str(error),
                }
                row["status"] = "skipped_shard_error"
                row["cell_count"] = 0
                row["failed_cell_count"] = None
                row["failures"] = [failure]
                if not args.skip_errors:
                    raise
            else:
                partial.replace(output)
                report["output"] = str(output)
                report["size_bytes"] = output.stat().st_size
                report["failures"] = []
                write_report(report, output.with_suffix(".report.json"))
                row["status"] = "written"
                row["cell_count"] = len(cells_and_paths)
                row["failed_cell_count"] = 0
                row["failures"] = []
                row["size_bytes"] = output.stat().st_size

        manifest["shards"].append(row)  # type: ignore[union-attr]
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(row), flush=True)

    print(json.dumps({"manifest": str(manifest_path), **manifest}, indent=2))


if __name__ == "__main__":
    main()
