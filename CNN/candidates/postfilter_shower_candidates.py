#!/usr/bin/env python3
"""Apply the Poisson-coherent shower-size cut to scan candidates."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-component-voxels", type=int, default=8)
    parser.add_argument("--sample", default="background_scan_candidate")
    args = parser.parse_args()
    if args.min_component_voxels < 1:
        raise ValueError("minimum component size must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with args.features.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = [row for row in reader if row["sample"] == args.sample]
    if not rows:
        raise ValueError(f"no rows found for sample {args.sample!r}")
    rows.sort(
        key=lambda row: (int(float(row["largest_component_voxels"])), float(row["q_hot"])),
        reverse=True,
    )
    accepted = [
        row for row in rows
        if int(float(row["largest_component_voxels"])) >= args.min_component_voxels
    ]
    rejected = [
        row for row in rows
        if int(float(row["largest_component_voxels"])) < args.min_component_voxels
    ]
    accepted_path = args.output_dir / f"candidates_size_ge{args.min_component_voxels}.csv"
    rejected_path = args.output_dir / f"candidates_size_lt{args.min_component_voxels}.csv"
    write_rows(accepted_path, accepted, fieldnames)
    write_rows(rejected_path, rejected, fieldnames)

    report = {
        "input": str(args.features),
        "sample": args.sample,
        "selection": f"largest_component_voxels >= {args.min_component_voxels}",
        "input_candidates": len(rows),
        "accepted_candidates": len(accepted),
        "rejected_candidates": len(rejected),
        "accepted_csv": str(accepted_path),
        "rejected_csv": str(rejected_path),
    }
    report_path = args.output_dir / f"postfilter_size_ge{args.min_component_voxels}.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
