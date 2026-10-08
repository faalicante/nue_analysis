#!/usr/bin/env python3
"""Collect sharded scan volumes into original event-group validation/test splits."""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import h5py
import numpy as np


SPLIT_NAMES = {0: "train", 1: "validation", 2: "test"}


def signal_split_map(reference_path: Path) -> dict[int, int]:
    with h5py.File(reference_path, "r") as source:
        signal = source["sample_type"][:] == 2
        event_ids = source["signal_event_id"][:][signal]
        splits = source["split"][:][signal]
    return {int(event): int(split) for event, split in zip(event_ids, splits)}


def create_like(output: h5py.File, source: h5py.File, path: str, count: int) -> h5py.Dataset:
    original = source[path]
    shape = (count, *original.shape[1:])
    kwargs: dict[str, object] = {}
    if path == "volumes_raw":
        kwargs.update(
            chunks=(1, 1, original.shape[2], original.shape[3]),
            compression="gzip", compression_opts=4, shuffle=True,
        )
    dataset = output.create_dataset(path, shape=shape, dtype=original.dtype, **kwargs)
    for key, value in original.attrs.items():
        dataset.attrs[key] = value
    return dataset


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--reference-hdf5", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--splits", nargs="+", choices=("validation", "test"), default=["validation", "test"])
    args = parser.parse_args()

    mapping = signal_split_map(args.reference_hdf5)
    requested_codes = {name: code for code, name in SPLIT_NAMES.items() if name in args.splits}
    input_paths = sorted(Path(path) for path in glob.glob(str(args.input_dir / "*.h5")))
    if not input_paths:
        raise ValueError(f"no HDF5 shards found in {args.input_dir}")

    locations: dict[str, list[tuple[Path, int]]] = {name: [] for name in requested_codes}
    for path in input_paths:
        with h5py.File(path, "r") as source:
            for index, event_id in enumerate(source["event_id"][:]):
                split = mapping.get(int(event_id))
                if split is not None and SPLIT_NAMES[split] in locations:
                    locations[SPLIT_NAMES[split]].append((path, index))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report: dict[str, object] = {"reference_hdf5": str(args.reference_hdf5), "splits": {}}
    with h5py.File(input_paths[0], "r") as template:
        dataset_paths: list[str] = []
        template.visititems(
            lambda name, obj: dataset_paths.append(name) if isinstance(obj, h5py.Dataset) else None
        )
        source_attrs = dict(template.attrs)

    for split_name, rows in locations.items():
        output_path = args.output_dir / f"signal_scan_{split_name}.h5"
        with h5py.File(output_path, "w") as output:
            for key, value in source_attrs.items():
                output.attrs[key] = value
            output.attrs["original_split"] = split_name
            with h5py.File(input_paths[0], "r") as template:
                outputs = {
                    path: create_like(output, template, path, len(rows))
                    for path in dataset_paths
                }
            current_path: Path | None = None
            current_source: h5py.File | None = None
            try:
                for output_index, (source_path, source_index) in enumerate(rows):
                    if source_path != current_path:
                        if current_source is not None:
                            current_source.close()
                        current_source = h5py.File(source_path, "r")
                        current_path = source_path
                    for path, dataset in outputs.items():
                        dataset[output_index] = current_source[path][source_index]
            finally:
                if current_source is not None:
                    current_source.close()
        report["splits"][split_name] = {  # type: ignore[index]
            "events": len(rows), "output": str(output_path), "size_bytes": output_path.stat().st_size,
        }
        print(json.dumps({"split": split_name, **report["splits"][split_name]}), flush=True)  # type: ignore[index]

    report_path = args.output_dir / "split_build_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(report_path), **report}, indent=2))


if __name__ == "__main__":
    main()
