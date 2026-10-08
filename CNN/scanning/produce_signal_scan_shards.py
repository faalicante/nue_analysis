#!/usr/bin/env python3
"""Produce restartable HDF5 shards of 57x200x200 signal scan volumes."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scanning.root_xypseg_to_scan_hdf5 import (
    CropSupportError,
    SignalTruth,
    read_truth_txt,
    validate_root_volume_support,
    write_hdf5,
)


LOOP_ENERGY_MIN_GEV = 80.0
LOOP_CROP_PLATE_MAX = 40
LOOP_THETA_MIN_MRAD = 5.0
LOOP_THETA_MAX_MRAD = 25.0
SCRIPT_VERSION = "2026-08-26-py39-1"


def theta_gt5_selection_failures(truth: SignalTruth) -> list[str]:
    """Apply only the lower-angle cut requested for scan-volume production."""
    theta_mrad = math.hypot(truth.ltx_mrad_corrected, truth.lty_mrad_corrected)
    return [] if theta_mrad > LOOP_THETA_MIN_MRAD else ["theta_not_gt_5_mrad"]


def training_selection_failures(truth: SignalTruth) -> list[str]:
    """Apply energy, plate and lower-theta cuts, with no upper theta cut."""
    failures = theta_gt5_selection_failures(truth)
    if not truth.energy_gev > LOOP_ENERGY_MIN_GEV:
        failures.append("energy_not_gt_80_gev")
    if not truth.plate - 3 < LOOP_CROP_PLATE_MAX:
        failures.append("crop_plate_not_lt_40")
    return failures


def loop_selection_failures(truth: SignalTruth) -> list[str]:
    """Return failed strict cuts, reproducing the selection in data_preparation/loop.py."""
    theta_mrad = math.hypot(truth.ltx_mrad_corrected, truth.lty_mrad_corrected)
    crop_plate = truth.plate - 3
    failures: list[str] = []
    if not truth.energy_gev > LOOP_ENERGY_MIN_GEV:
        failures.append("energy_not_gt_80_gev")
    if not crop_plate < LOOP_CROP_PLATE_MAX:
        failures.append("crop_plate_not_lt_40")
    if not theta_mrad > LOOP_THETA_MIN_MRAD:
        failures.append("theta_not_gt_5_mrad")
    if not theta_mrad < LOOP_THETA_MAX_MRAD:
        failures.append("theta_not_lt_25_mrad")
    return failures


def read_event_list(path: Path) -> set[int]:
    events: set[int] = set()
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        value = line.split("#", 1)[0].strip()
        if not value:
            continue
        try:
            events.add(int(value))
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: invalid event id {value!r}") from error
    return events


def render_root_path(template: str, event_id: int) -> str:
    try:
        return template.format(event=event_id, event_plus_one=event_id + 1)
    except KeyError as error:
        raise ValueError(
            "root template may only use {event} and {event_plus_one}"
        ) from error


def chunked(values: list[int], size: int) -> list[list[int]]:
    if size <= 0:
        raise ValueError("shard size must be positive")
    return [values[first : first + size] for first in range(0, len(values), size)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth-txt", type=Path, required=True)
    parser.add_argument(
        "--root-template",
        required=True,
        help="For example /eos/.../b000021.0.0.{event_plus_one}.trk.root",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-size", type=int, default=100)
    parser.add_argument("--event-min", type=int)
    parser.add_argument("--event-max", type=int)
    parser.add_argument("--events-file", type=Path)
    parser.add_argument("--max-events", type=int)
    parser.add_argument(
        "--selection", choices=("training", "theta_gt5", "loop", "none"), default="training",
        help="Default: energy >80, plate-3 <40, theta >5, with no upper theta cut",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    truth = read_truth_txt(args.truth_txt)
    selected = sorted(truth)
    rejection_counts: dict[str, int] = {}
    if args.selection != "none":
        accepted: list[int] = []
        for event in selected:
            if args.selection == "training":
                failures = training_selection_failures(truth[event])
            elif args.selection == "theta_gt5":
                failures = theta_gt5_selection_failures(truth[event])
            else:
                failures = loop_selection_failures(truth[event])
            if failures:
                for reason in failures:
                    rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
            else:
                accepted.append(event)
        selected = accepted
    if args.event_min is not None:
        selected = [event for event in selected if event >= args.event_min]
    if args.event_max is not None:
        selected = [event for event in selected if event <= args.event_max]
    if args.events_file is not None:
        requested = read_event_list(args.events_file)
        missing = sorted(requested - set(truth))
        if missing:
            raise ValueError(f"events absent from truth table: {missing[:10]}")
        selected = [event for event in selected if event in requested]
    if args.max_events is not None:
        if args.max_events <= 0:
            raise ValueError("max-events must be positive")
        selected = selected[: args.max_events]
    if not selected:
        raise ValueError("event selection is empty")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {
        "producer_version": SCRIPT_VERSION,
        "truth_txt": str(args.truth_txt),
        "root_template": args.root_template,
        "selection": {
            "name": args.selection,
            "strict_cuts": None if args.selection == "none" else ({
                "energy_gev": "> 80",
                "crop_plate": "truth_plate - 3 < 40",
                "theta_mrad": "> 5 using hypot(lt_corrected_x, lt_corrected_y)",
                "theta_max_mrad": None,
                "slope_correction": "lt_corrected = lt_raw * 1315 / 1350",
            } if args.selection == "training" else {
                "theta_mrad": "> 5 using hypot(lt_corrected_x, lt_corrected_y)",
                "slope_correction": "lt_corrected = lt_raw * 1315 / 1350",
            } if args.selection == "theta_gt5" else {
                "energy_gev": "> 80",
                "crop_plate": "truth_plate - 3 < 40",
                "theta_mrad": "5 < hypot(lt_corrected_x, lt_corrected_y) < 25",
                "slope_correction": "lt_corrected = lt_raw * 1315 / 1350",
            }),
            "input_truth_events": len(truth),
            "rejection_counts_nonexclusive": rejection_counts,
        },
        "shard_size": args.shard_size,
        "selected_events": len(selected),
        "shards": [],
    }
    for shard_index, event_ids in enumerate(chunked(selected, args.shard_size)):
        output = args.output_dir / (
            f"signal_scan_{shard_index:04d}_evt{event_ids[0]:06d}-{event_ids[-1]:06d}.h5"
        )
        root_paths = [render_root_path(args.root_template, event) for event in event_ids]
        row: dict[str, object] = {
            "index": shard_index,
            "output": str(output),
            "event_first": event_ids[0],
            "event_last": event_ids[-1],
            "event_count": len(event_ids),
        }
        if output.exists() and not args.overwrite:
            row["status"] = "skipped_existing"
        elif args.dry_run:
            row["status"] = "dry_run"
            row["first_root"] = root_paths[0]
        else:
            usable_paths: list[str] = []
            skipped_support: list[dict[str, object]] = []
            for event_id, root_path in zip(event_ids, root_paths):
                try:
                    validate_root_volume_support(root_path, truth[event_id])
                except CropSupportError as error:
                    skipped_support.append({
                        "event_id": event_id,
                        "source_file": root_path,
                        "reason": str(error),
                    })
                    print(json.dumps({
                        "status": "skipped_insufficient_support",
                        "event_id": event_id,
                        "reason": str(error),
                    }), flush=True)
                else:
                    usable_paths.append(root_path)
            row["requested_event_count"] = len(event_ids)
            row["skipped_insufficient_support"] = len(skipped_support)
            if not usable_paths:
                row["event_count"] = 0
                row["status"] = "skipped_no_usable_events"
                manifest["shards"].append(row)  # type: ignore[union-attr]
                manifest_path = args.output_dir / "manifest.json"
                manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
                print(json.dumps(row), flush=True)
                continue
            partial = output.with_suffix(".partial.h5")
            if partial.exists():
                partial.unlink()
            report = write_hdf5(usable_paths, truth, partial, reference_samples=None)
            partial.replace(output)
            report["output"] = str(output)
            report["size_bytes"] = output.stat().st_size
            report["skipped_insufficient_support"] = skipped_support
            output.with_suffix(".report.json").write_text(
                json.dumps(report, indent=2) + "\n", encoding="utf-8"
            )
            row["status"] = "written"
            row["event_count"] = len(usable_paths)
            row["size_bytes"] = report["size_bytes"]
        manifest["shards"].append(row)  # type: ignore[union-attr]
        manifest_path = args.output_dir / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(row), flush=True)

    print(json.dumps({"manifest": str(args.output_dir / "manifest.json"), **manifest}, indent=2))


if __name__ == "__main__":
    main()
