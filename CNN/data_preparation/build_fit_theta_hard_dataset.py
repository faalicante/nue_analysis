#!/usr/bin/env python3
"""Build a training HDF5 with hard crops selected by coherent Poisson-hot fit."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path

import h5py
import numpy as np
import uproot


SCALAR_DATASETS = (
    "presence",
    "sample_type",
    "slope_xy",
    "signal_event_id",
    "tile_id",
    "crop_id_in_tile",
    "negative_type",
    "split",
    "source_entry",
)


def selected_entries(path: Path, theta_max_mrad: float) -> set[int]:
    selected: set[int] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["fit_quality"] != "coherent" or not row["fit_theta_mrad"].strip():
                continue
            if float(row["fit_theta_mrad"]) < theta_max_mrad:
                selected.add(int(row["root_entry"]))
    if not selected:
        raise ValueError(f"no coherent entries below {theta_max_mrad} mrad in {path}")
    return selected


def root_scalars(path: Path, entries: np.ndarray) -> dict[str, np.ndarray]:
    with uproot.open(path) as root_file:
        tree = root_file["samples"]
        arrays = tree.arrays(
            ["background_mu", "presence", "sample_type", "slope_x", "slope_y", "signal_event_id", "cell_id", "tag_cell_id", "zScale"],
            library="np",
        )
    result = {name: np.asarray(values)[entries] for name, values in arrays.items()}
    if not np.all(result["sample_type"] == 1) or not np.all(result["presence"] == 1):
        raise ValueError(f"{path} selected rows are not all hard crops")
    if not np.all(result["signal_event_id"] == -1):
        raise ValueError(f"{path} hard rows have a signal event id")
    if not np.allclose(result["slope_x"], 0.0) or not np.allclose(result["slope_y"], 0.0):
        raise ValueError(f"{path} hard rows have non-zero slope labels")
    return result


def copy_regions(
    source: h5py.Dataset,
    keep: np.ndarray,
    target: h5py.Dataset,
    offset: int,
    chunk_size: int,
) -> int:
    for start in range(0, len(keep), chunk_size):
        stop = min(start + chunk_size, len(keep))
        local_keep = keep[start:stop]
        if not np.any(local_keep):
            continue
        values = np.asarray(source[start:stop])[local_keep]
        target[offset : offset + len(values)] = values
        offset += len(values)
    return offset


def append_root_regions(
    path: Path,
    entries: np.ndarray,
    target: h5py.Dataset,
    offset: int,
    chunk_size: int,
) -> int:
    entry_set = set(int(value) for value in entries)
    with uproot.open(path) as root_file:
        tree = root_file["samples"]
        for start in range(0, int(tree.num_entries), chunk_size):
            stop = min(start + chunk_size, int(tree.num_entries))
            local_entries = np.arange(start, stop, dtype=np.int64)
            local_keep = np.fromiter(
                (int(entry) in entry_set for entry in local_entries), dtype=bool, count=len(local_entries)
            )
            if not np.any(local_keep):
                continue
            values = np.asarray(tree["counts"].array(entry_start=start, entry_stop=stop, library="np"))
            target[offset : offset + int(np.count_nonzero(local_keep))] = values[local_keep]
            offset += int(np.count_nonzero(local_keep))
    return offset


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-hdf5", type=Path, required=True)
    parser.add_argument("--bkg-root", type=Path, required=True)
    parser.add_argument("--bkg15-root", type=Path, required=True)
    parser.add_argument("--bkg-angles", type=Path, required=True)
    parser.add_argument("--bkg15-angles", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--theta-max-mrad", type=float, default=5.0)
    parser.add_argument("--exclude-cell-id", type=int, action="append", default=[58])
    parser.add_argument("--chunk-size", type=int, default=16)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.theta_max_mrad <= 0.0 or args.chunk_size <= 0:
        raise ValueError("theta-max-mrad and chunk-size must be positive")
    output = args.output.expanduser().resolve()
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"{output} already exists; use --overwrite")
    excluded_cells = set(int(value) for value in args.exclude_cell_id)
    selected_bkg = selected_entries(args.bkg_angles.expanduser(), args.theta_max_mrad)
    selected_bkg15 = selected_entries(args.bkg15_angles.expanduser(), args.theta_max_mrad)

    base_path = args.base_hdf5.expanduser().resolve()
    bkg_path = str(args.bkg_root.expanduser().resolve())
    bkg15_path = args.bkg15_root.expanduser().resolve()
    with h5py.File(base_path, "r") as base:
        source_files = base["source_file"].asstr()[:]
        sample_type = np.asarray(base["sample_type"][:], dtype=np.int8)
        source_entry = np.asarray(base["source_entry"][:], dtype=np.int64)
        tile_id = np.asarray(base["tile_id"][:], dtype=np.int64)
        is_original_hard = sample_type == 1
        if not np.all(source_files[is_original_hard] == bkg_path):
            raise ValueError("base HDF5 hard rows do not all originate from --bkg-root")
        keep = ~is_original_hard | np.isin(source_entry, list(selected_bkg))
        keep &= ~((sample_type == 1) & np.isin(tile_id, list(excluded_cells)))
        kept_hard = int(np.count_nonzero(keep & is_original_hard))

        split_by_tile = {
            int(cell): int(code)
            for cell, code in zip(tile_id[is_original_hard], np.asarray(base["split"][:])[is_original_hard], strict=True)
        }

        bkg15_entries = np.asarray(sorted(selected_bkg15), dtype=np.int64)
        bkg15_scalars = root_scalars(bkg15_path, bkg15_entries)
        append_keep = ~np.isin(bkg15_scalars["cell_id"], list(excluded_cells))
        bkg15_entries = bkg15_entries[append_keep]
        bkg15_scalars = {name: values[append_keep] for name, values in bkg15_scalars.items()}
        missing_split_cells = sorted(set(int(cell) for cell in bkg15_scalars["cell_id"]) - set(split_by_tile))
        if missing_split_cells:
            raise ValueError(f"no preserved split assignment for bkg_15 cells {missing_split_cells}")
        append_split = np.asarray(
            [split_by_tile[int(cell)] for cell in bkg15_scalars["cell_id"]], dtype=np.int8
        )
        total = int(np.count_nonzero(keep)) + len(bkg15_entries)

        output.parent.mkdir(parents=True, exist_ok=True)
        temporary_handle = tempfile.NamedTemporaryFile(
            prefix=f".{output.name}.", suffix=".tmp", dir=output.parent, delete=False
        )
        temporary = Path(temporary_handle.name)
        temporary_handle.close()
        try:
            with h5py.File(temporary, "w") as target:
                for key, value in base.attrs.items():
                    target.attrs[key] = value
                target.attrs["hard_selection"] = (
                    f"Poisson-hot fit_quality=coherent and fit_theta_mrad<{args.theta_max_mrad}"
                )
                target.attrs["hard_selection_excluded_cell_ids_json"] = json.dumps(sorted(excluded_cells))
                target.attrs["hard_selection_added_source"] = str(bkg15_path)
                regions = target.create_dataset(
                    "regions_raw", shape=(total, 57, 32, 32), dtype=base["regions_raw"].dtype,
                    chunks=(min(args.chunk_size, total), 57, 32, 32), compression="gzip", compression_opts=4,
                    shuffle=True,
                )
                offset = copy_regions(base["regions_raw"], keep, regions, 0, args.chunk_size)
                offset = append_root_regions(bkg15_path, bkg15_entries, regions, offset, args.chunk_size)
                if offset != total:
                    raise RuntimeError(f"wrote {offset} region rows, expected {total}")

                for name in SCALAR_DATASETS:
                    base_values = np.asarray(base[name][:])[keep]
                    if name == "presence":
                        appended = bkg15_scalars["presence"].astype(base_values.dtype)
                    elif name == "sample_type":
                        appended = bkg15_scalars["sample_type"].astype(base_values.dtype)
                    elif name == "slope_xy":
                        appended = np.zeros((len(bkg15_entries), 2), dtype=base_values.dtype)
                    elif name == "signal_event_id":
                        appended = bkg15_scalars["signal_event_id"].astype(base_values.dtype)
                    elif name == "tile_id":
                        appended = bkg15_scalars["cell_id"].astype(base_values.dtype)
                    elif name == "crop_id_in_tile":
                        appended = bkg15_scalars["tag_cell_id"].astype(base_values.dtype)
                    elif name == "negative_type":
                        appended = np.ones(len(bkg15_entries), dtype=base_values.dtype)
                    elif name == "split":
                        appended = append_split.astype(base_values.dtype)
                    else:
                        appended = bkg15_entries.astype(base_values.dtype)
                    target.create_dataset(name, data=np.concatenate((base_values, appended)))

                base_background = np.asarray(base["background_mu"][:])
                if base_background.ndim != 2:
                    raise ValueError("base HDF5 must store one background_mu value per event and z")
                appended_background = np.repeat(
                    bkg15_scalars["background_mu"].astype(base_background.dtype)[:, None], 57, axis=1
                )
                target.create_dataset("background_mu", data=np.concatenate((base_background[keep], appended_background)))

                string_dtype = h5py.string_dtype(encoding="utf-8")
                appended_sources = np.full(len(bkg15_entries), str(bkg15_path), dtype=object)
                target.create_dataset(
                    "source_file", data=np.concatenate((source_files[keep], appended_sources)), dtype=string_dtype
                )
                metadata = target.create_group("metadata")
                if "metadata/zScale" in base:
                    zscale = np.concatenate((
                        np.asarray(base["metadata/zScale"][:])[keep],
                        bkg15_scalars["zScale"].astype(np.asarray(base["metadata/zScale"][:]).dtype),
                    ))
                    metadata.create_dataset("zScale", data=zscale)
                split_group = target.create_group("split_indices")
                combined_split = np.asarray(target["split"][:], dtype=np.int8)
                for code, name in enumerate(("train", "validation", "test")):
                    split_group.create_dataset(name, data=np.flatnonzero(combined_split == code).astype(np.int64))
            os.replace(temporary, output)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise

    summary = {
        "output": str(output),
        "theta_max_mrad": args.theta_max_mrad,
        "fit_quality": "coherent",
        "excluded_cell_ids": sorted(excluded_cells),
        "bkg_selected_from_angles": len(selected_bkg),
        "bkg_selected_after_existing_cell_filter": kept_hard,
        "bkg15_selected_from_angles": len(selected_bkg15),
        "bkg15_added": len(bkg15_entries),
        "hard_total": kept_hard + len(bkg15_entries),
        "total_examples": total,
    }
    summary_path = output.with_suffix(".selection_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
