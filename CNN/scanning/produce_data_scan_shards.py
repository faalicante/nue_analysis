#!/usr/bin/env python3
"""Produce restartable 57x200x200 HDF5 shards for real-data CNN scans.

The input layout is configurable because the directory above ``cell_X0_Y0``
can differ between data productions.  The ROOT template may use ``{rwb}``
(three-digit run-wall-brick token), ``{brick}`` (for example ``b000021``),
``{x}``, ``{y}``, and ``{cell_id}``.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from scanning.produce_background_scan_shards import chunked
from scanning.root_background_xypseg_to_scan_hdf5 import (
    BackgroundCell,
    CropSupportError,
    validate_background_root,
    validate_partial_background_root,
    write_background_hdf5,
    write_partial_background_hdf5,
    write_report,
)


SCRIPT_VERSION = "2026-09-08-1"


DATA_CELL_SIZE_UM = 10_000.0


@dataclass(frozen=True)
class DataCell(BackgroundCell):
    """Data cell whose directory coordinate is its center in units of 10 mm."""

    @property
    def x_bounds_um(self) -> tuple[float, float]:
        center = self.x * DATA_CELL_SIZE_UM
        return center - DATA_CELL_SIZE_UM / 2, center + DATA_CELL_SIZE_UM / 2

    @property
    def y_bounds_um(self) -> tuple[float, float]:
        center = self.y * DATA_CELL_SIZE_UM
        return center - DATA_CELL_SIZE_UM / 2, center + DATA_CELL_SIZE_UM / 2


def selected_data_cells(x_min: int, x_max: int, y_min: int, y_max: int) -> list[DataCell]:
    if not (1 <= x_min <= x_max <= 18 and 1 <= y_min <= y_max <= 18):
        raise ValueError("x/y ranges must satisfy 1 <= min <= max <= 18")
    return [
        DataCell(x=x, y=y)
        for y in range(y_min, y_max + 1)
        for x in range(x_min, x_max + 1)
    ]


def cell_from_id(cell_id: int) -> DataCell:
    if not 0 <= cell_id < 18 * 18:
        raise ValueError(f"cell_id must be in [0, 323], got {cell_id}")
    return DataCell(x=cell_id % 18 + 1, y=cell_id // 18 + 1)


def read_cell_list(path: Path) -> list[DataCell]:
    """Read selected cells while preserving order and removing duplicates.

    Accepted records are ``cell_id``, ``x y``, ``x,y``, ``cell_id x y``, or a
    path/name containing ``cell_X0_Y0``.  Text after ``#`` is a comment.
    """
    cells: list[DataCell] = []
    seen: set[int] = set()
    cell_name = re.compile(r"cell_([1-9]|1[0-8])0_([1-9]|1[0-8])0(?:\D|$)")
    for line_number, original in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        match = cell_name.search(line)
        try:
            if match:
                cell = DataCell(x=int(match.group(1)), y=int(match.group(2)))
            else:
                fields = line.replace(",", " ").split()
                if fields and all(field.lower() in {"cell_id", "cell", "id", "x", "y"} for field in fields):
                    continue
                values = [int(field) for field in fields]
                if len(values) == 1:
                    cell = cell_from_id(values[0])
                elif len(values) == 2:
                    cell = DataCell(x=values[0], y=values[1])
                elif len(values) == 3:
                    expected = DataCell(x=values[1], y=values[2])
                    if values[0] != expected.cell_id:
                        raise ValueError(
                            f"cell_id {values[0]} disagrees with x={values[1]}, y={values[2]}"
                        )
                    cell = expected
                else:
                    raise ValueError("expected 1, 2, or 3 integer columns")
        except (TypeError, ValueError) as error:
            raise ValueError(f"{path}:{line_number}: cannot parse {original!r}: {error}") from error
        if cell.cell_id not in seen:
            cells.append(cell)
            seen.add(cell.cell_id)
    if not cells:
        raise ValueError(f"cell list {path} contains no usable cells")
    return cells


def normalize_rwb(value: str) -> tuple[int, str, str]:
    """Return numeric id, zero-padded token and FEDRA brick directory name."""
    token = value.strip()
    if token.startswith("b000"):
        token = token[4:]
    if not token.isdigit():
        raise ValueError(f"rwb must be numeric or have form b000NNN, got {value!r}")
    numeric = int(token)
    if numeric < 0:
        raise ValueError("rwb must be non-negative")
    padded = f"{numeric:03d}"
    return numeric, padded, f"b000{padded}"


def render_data_root_path(
    template: str,
    cell: DataCell,
    *,
    rwb: str,
    brick: str,
) -> str:
    try:
        return template.format(
            rwb=rwb,
            brick=brick,
            x=cell.x,
            y=cell.y,
            cell_id=cell.cell_id,
        )
    except KeyError as error:
        raise ValueError(
            "root template may only use {rwb}, {brick}, {x}, {y}, and {cell_id}"
        ) from error


def add_data_provenance(
    path: Path,
    *,
    rwb_id: int,
    rwb: str,
    brick: str,
    sample_kind: str = "data_full_cell",
) -> None:
    """Label a data shard and make its event ids globally unique."""
    with h5py.File(path, "r+") as hdf5:
        count = int(hdf5["event_id"].shape[0])
        cell_ids = np.asarray(hdf5["cell_id"][:], dtype=np.int64)
        # There are at most 324 cells, so a base of 1000 is collision-free.
        hdf5["event_id"][:] = np.int64(rwb_id) * 1000 + cell_ids
        string_dtype = h5py.string_dtype(encoding="utf-8")
        hdf5.create_dataset(
            "run_wall_brick", data=np.full(count, rwb_id, dtype=np.int64)
        )
        hdf5.create_dataset(
            "rwb_token", data=np.asarray([rwb] * count, dtype=object), dtype=string_dtype
        )
        hdf5.create_dataset(
            "brick_name", data=np.asarray([brick] * count, dtype=object), dtype=string_dtype
        )
        hdf5.attrs["sample_kind"] = sample_kind
        hdf5.attrs["event_id_formula"] = "run_wall_brick * 1000 + cell_id"
        hdf5.attrs["run_wall_brick"] = rwb_id
        hdf5.attrs["brick_name"] = brick
        hdf5.attrs["data_cell_center_formula_um"] = "x*10000, y*10000"
        hdf5.attrs["data_cell_bounds_formula_um"] = "(coordinate-0.5)*10000 .. (coordinate+0.5)*10000"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert per-cell data ROOT files into full-cell CNN scan shards."
    )
    parser.add_argument(
        "--root-template",
        required=True,
        help=(
            "ROOT path template; e.g. '/eos/.../cell_{x}0_{y}0/{brick}/"
            "{brick}.0.0.0.trk.root'"
        ),
    )
    parser.add_argument(
        "--rwb",
        required=True,
        help="Run-wall-brick number, e.g. 21, 021, or b000021",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-size", type=int, default=25)
    parser.add_argument("--x-min", type=int, default=1)
    parser.add_argument("--x-max", type=int, default=18)
    parser.add_argument("--y-min", type=int, default=1)
    parser.add_argument("--y-max", type=int, default=18)
    parser.add_argument(
        "--cell-list",
        type=Path,
        help="Use only cells listed in this text file; overrides the x/y rectangle",
    )
    parser.add_argument("--max-cells", type=int)
    parser.add_argument(
        "--exclude-cell-id",
        action="append",
        default=[],
        type=int,
        help="Cell id to exclude; repeat the option for multiple cells",
    )
    parser.add_argument(
        "--exclude-cell-list",
        type=Path,
        help="Text file of cells to exclude, using the same formats as --cell-list",
    )
    parser.add_argument(
        "--skip-errors",
        action="store_true",
        help="Record missing/invalid ROOT files and continue; default is fail-fast",
    )
    parser.add_argument(
        "--allow-partial-cells",
        action="store_true",
        help=(
            "Export a clipped cell as its own scan HDF5 when it cannot support "
            "the usual 200x200 (1 cm) volume but still contains at least 20x20 bins."
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    rwb_id, rwb, brick = normalize_rwb(args.rwb)
    cells = (
        read_cell_list(args.cell_list)
        if args.cell_list is not None
        else selected_data_cells(args.x_min, args.x_max, args.y_min, args.y_max)
    )
    excluded_cell_ids = set(args.exclude_cell_id)
    if args.exclude_cell_list is not None:
        excluded_cell_ids.update(
            cell.cell_id for cell in read_cell_list(args.exclude_cell_list)
        )
    excluded_cell_ids = sorted(excluded_cell_ids)
    if any(cell_id < 0 or cell_id >= 18 * 18 for cell_id in excluded_cell_ids):
        raise ValueError("exclude-cell-id values must be in [0, 323]")
    cells = [cell for cell in cells if cell.cell_id not in excluded_cell_ids]
    if args.max_cells is not None:
        if args.max_cells <= 0:
            raise ValueError("max-cells must be positive")
        cells = cells[: args.max_cells]
    if not cells:
        raise ValueError("cell selection is empty")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / f"manifest_{brick}.json"
    manifest: dict[str, object] = {
        "producer_version": SCRIPT_VERSION,
        "sample_kind": "data_full_cell",
        "run_wall_brick": rwb_id,
        "rwb_token": rwb,
        "brick_name": brick,
        "root_template": args.root_template,
        "cell_list": None if args.cell_list is None else str(args.cell_list.resolve()),
        "exclude_cell_list": (
            None if args.exclude_cell_list is None else str(args.exclude_cell_list.resolve())
        ),
        "cell_ranges": {"x": [args.x_min, args.x_max], "y": [args.y_min, args.y_max]},
        "excluded_cell_ids": excluded_cell_ids,
        "selected_cells": len(cells),
        "shard_size": args.shard_size,
        "allow_partial_cells": args.allow_partial_cells,
        "shards": [],
    }

    for shard_index, shard_cells in enumerate(chunked(cells, args.shard_size)):
        first, last = shard_cells[0], shard_cells[-1]
        output = args.output_dir / (
            f"data_scan_{brick}_{shard_index:04d}_"
            f"cell{first.cell_id:03d}-{last.cell_id:03d}.h5"
        )
        row: dict[str, object] = {
            "index": shard_index,
            "output": str(output),
            "cell_id_first": first.cell_id,
            "cell_id_last": last.cell_id,
            "requested_cell_count": len(shard_cells),
        }
        if output.exists() and not args.overwrite and not args.allow_partial_cells:
            row["status"] = "skipped_existing"
            row["cell_count"] = len(shard_cells)
        elif args.dry_run:
            row["status"] = "dry_run"
            row["cell_count"] = len(shard_cells)
            row["first_root"] = render_data_root_path(
                args.root_template, first, rwb=rwb, brick=brick
            )
        else:
            usable: list[tuple[BackgroundCell, str]] = []
            partial_usable: list[tuple[BackgroundCell, str]] = []
            failures: list[dict[str, object]] = []
            for cell in shard_cells:
                root_path = render_data_root_path(
                    args.root_template, cell, rwb=rwb, brick=brick
                )
                try:
                    validate_background_root(root_path, cell)
                except CropSupportError as error:
                    if args.allow_partial_cells:
                        try:
                            validate_partial_background_root(root_path, cell)
                        except Exception as partial_error:
                            error = partial_error
                        else:
                            partial_usable.append((cell, root_path))
                            print(json.dumps({
                                "status": "partial_cell",
                                "cell_id": cell.cell_id,
                                "cell_x": cell.x,
                                "cell_y": cell.y,
                                "source_file": root_path,
                                "reason": str(error),
                            }), flush=True)
                            continue
                    failure = {
                        "cell_id": cell.cell_id,
                        "cell_x": cell.x,
                        "cell_y": cell.y,
                        "source_file": root_path,
                        "error_type": type(error).__name__,
                        "reason": str(error),
                    }
                    failures.append(failure)
                    print(json.dumps({"status": "cell_error", **failure}), flush=True)
                    if not args.skip_errors:
                        raise error
                except Exception as error:
                    failure = {
                        "cell_id": cell.cell_id,
                        "cell_x": cell.x,
                        "cell_y": cell.y,
                        "source_file": root_path,
                        "error_type": type(error).__name__,
                        "reason": str(error),
                    }
                    failures.append(failure)
                    print(json.dumps({"status": "cell_error", **failure}), flush=True)
                    if not args.skip_errors:
                        raise
                else:
                    usable.append((cell, root_path))

            row["failed_cell_count"] = len(failures)
            row["failures"] = failures
            partial_reports: list[dict[str, object]] = []
            for cell, root_path in partial_usable:
                partial_index = 5000 + cell.cell_id
                partial_output = args.output_dir / (
                    f"data_scan_{brick}_{partial_index:04d}_"
                    f"cell{cell.cell_id:03d}-{cell.cell_id:03d}.h5"
                )
                if partial_output.exists() and not args.overwrite:
                    partial_report: dict[str, object] = {
                        "status": "skipped_existing",
                        "output": str(partial_output),
                        "cell_id": cell.cell_id,
                        "cell_x": cell.x,
                        "cell_y": cell.y,
                        "partial_cell": True,
                    }
                else:
                    partial = partial_output.with_suffix(".partial.h5")
                    if partial.exists():
                        partial.unlink()
                    partial_report = write_partial_background_hdf5(cell, root_path, partial)
                    add_data_provenance(
                        partial,
                        rwb_id=rwb_id,
                        rwb=rwb,
                        brick=brick,
                        sample_kind="data_partial_cell",
                    )
                    partial.replace(partial_output)
                    partial_report.update({
                        "status": "written",
                        "output": str(partial_output),
                        "size_bytes": partial_output.stat().st_size,
                    })
                    write_report(partial_report, partial_output.with_suffix(".report.json"))
                partial_reports.append(partial_report)
            row["partial_cell_count"] = len(partial_reports)
            row["partial_cells"] = partial_reports
            if not usable and not partial_usable:
                row["status"] = "skipped_no_usable_cells"
                row["cell_count"] = 0
            elif usable:
                if output.exists() and not args.overwrite:
                    row["status"] = "skipped_existing_with_partial_cells"
                    row["cell_count"] = len(usable)
                    row["size_bytes"] = output.stat().st_size
                    manifest["shards"].append(row)  # type: ignore[union-attr]
                    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
                    print(json.dumps(row), flush=True)
                    continue
                partial = output.with_suffix(".partial.h5")
                if partial.exists():
                    partial.unlink()
                report = write_background_hdf5(usable, partial)
                add_data_provenance(
                    partial, rwb_id=rwb_id, rwb=rwb, brick=brick
                )
                partial.replace(output)
                report.update(
                    {
                        "sample_kind": "data_full_cell",
                        "run_wall_brick": rwb_id,
                        "rwb_token": rwb,
                        "brick_name": brick,
                        "output": str(output),
                        "size_bytes": output.stat().st_size,
                        "failures": failures,
                    }
                )
                write_report(report, output.with_suffix(".report.json"))
                row["status"] = "written"
                row["cell_count"] = len(usable)
                row["size_bytes"] = output.stat().st_size
            else:
                row["status"] = "written_partial_only"
                row["cell_count"] = 0

        manifest["shards"].append(row)  # type: ignore[union-attr]
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(row), flush=True)

    print(json.dumps({"manifest": str(manifest_path), **manifest}, indent=2))


if __name__ == "__main__":
    main()
