#!/usr/bin/env python3
"""Align HDF5 train/validation/test groups to an existing manifest.

Rows are matched by the same leakage-prevention keys used by data_preparation/root_to_hdf5.py:
signal_event_id for signal and tile_id for background samples.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


SPLITS = ("train", "validation", "test")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--statistics", type=Path, required=True)
    args = parser.parse_args()

    reference: dict[tuple[str, int], str] = {}
    with args.reference_manifest.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if int(row["sample_type"]) == 2:
                key = ("signal_event_id", int(row["signal_event_id"]))
            else:
                key = ("tile_id", int(row["tile_id"]))
            previous = reference.setdefault(key, row["split"])
            if previous != row["split"]:
                raise ValueError(f"reference split leakage for {key}: {previous} vs {row['split']}")

    with h5py.File(args.hdf5, "r+") as hdf5:
        sample_type = np.asarray(hdf5["sample_type"])
        signal_id = np.asarray(hdf5["signal_event_id"])
        tile_id = np.asarray(hdf5["tile_id"])
        split_names = [
            reference[("signal_event_id", int(event_id))]
            if int(kind) == 2 else reference[("tile_id", int(cell_id))]
            for kind, event_id, cell_id in zip(sample_type, signal_id, tile_id, strict=True)
        ]
        split = np.asarray([SPLITS.index(name) for name in split_names], dtype=np.int8)
        hdf5["split"][:] = split
        group = hdf5["split_indices"]
        for code, name in enumerate(SPLITS):
            del group[name]
            group.create_dataset(name, data=np.flatnonzero(split == code).astype(np.int64))

    manifest_rows: list[dict[str, str]] = []
    with args.manifest.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fieldnames = reader.fieldnames
        assert fieldnames is not None
        for index, row in enumerate(reader):
            row["split"] = split_names[index]
            manifest_rows.append(row)
    temporary_manifest = args.manifest.with_suffix(args.manifest.suffix + ".tmp")
    with temporary_manifest.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)
    temporary_manifest.replace(args.manifest)

    statistics = json.loads(args.statistics.read_text(encoding="utf-8"))
    statistics["split"] = {}
    for code, name in enumerate(SPLITS):
        mask = split == code
        statistics["split"][name] = {
            "total": int(mask.sum()),
            "poisson": int(np.count_nonzero(mask & (sample_type == 0))),
            "hard": int(np.count_nonzero(mask & (sample_type == 1))),
            "signal": int(np.count_nonzero(mask & (sample_type == 2))),
        }
    args.statistics.write_text(json.dumps(statistics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"hdf5": str(args.hdf5), "split": statistics["split"]}, indent=2))


if __name__ == "__main__":
    main()
