#!/usr/bin/env python3
"""Create smoothed layer GIFs for ROOT crops selected from an angle CSV."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import uproot

from visualization.make_background_candidate_gifs import root_th2_smooth
from visualization.make_original_crop_gifs import render_frames
from visualization.make_root_crop_xzyz_projections import number, selected_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--angles", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--angle-source", choices=("cnn", "fit"), default="cnn")
    parser.add_argument("--min-angle-mrad", type=float, default=10.0)
    parser.add_argument("--score-min", type=float, default=0.0)
    parser.add_argument("--coherent-fit-only", action="store_true")
    parser.add_argument("--count", type=int, default=20)
    parser.add_argument("--frame-duration-ms", type=int, default=160)
    args = parser.parse_args()
    if args.min_angle_mrad < 0.0:
        parser.error("--min-angle-mrad must be non-negative")
    if not 0.0 <= args.score_min <= 1.0:
        parser.error("--score-min must be in [0, 1]")
    if args.count <= 0:
        parser.error("--count must be positive")

    angle_field = f"{args.angle_source}_theta_mrad"
    rows = selected_rows(
        args.angles, angle_field=angle_field, min_angle=args.min_angle_mrad,
        score_min=args.score_min, coherent_fit_only=args.coherent_fit_only,
        limit=args.count,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    with uproot.open(args.input_root) as root_file:
        tree = root_file["samples"]
        for rank, row in enumerate(rows, start=1):
            entry = int(row["root_entry"])
            arrays = tree.arrays(
                ["counts", "background_mu"], entry_start=entry, entry_stop=entry + 1, library="np"
            )
            raw = np.asarray(arrays["counts"][0], dtype=np.float32)
            title = (
                f"ROOT entry {entry} · cell={row['cell_id']} · tag={row['tag_cell_id']}\n"
                f"score={number(row, 'cnn_presence_score'):.3f} · "
                f"theta CNN={number(row, 'cnn_theta_mrad'):.1f} mrad · "
                f"theta fit={number(row, 'fit_theta_mrad'):.1f} mrad · {row['fit_quality']}"
            )
            output = args.output_dir / (
                f"rank_{rank:03d}_entry_{entry:05d}_theta_{number(row, angle_field):06.2f}.gif"
            )
            frames = render_frames(root_th2_smooth(raw), title)
            durations = [1000] + [args.frame_duration_ms] * 55 + [1000]
            frames[0].save(
                output, save_all=True, append_images=frames[1:], duration=durations,
                loop=0, optimize=False, disposal=2,
            )
            manifest_rows.append({
                "rank": rank,
                "angle_source": args.angle_source,
                "selected_angle_mrad": number(row, angle_field),
                "root_entry": entry,
                "cell_id": int(row["cell_id"]),
                "tag_cell_id": int(row["tag_cell_id"]),
                "cnn_presence_score": number(row, "cnn_presence_score"),
                "cnn_theta_mrad": number(row, "cnn_theta_mrad"),
                "fit_theta_mrad": number(row, "fit_theta_mrad"),
                "fit_quality": row["fit_quality"],
                "gif": str(output.resolve()),
            })
            print(json.dumps({"rank": rank, "root_entry": entry, "gif": str(output)}), flush=True)
    with (args.output_dir / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = list(manifest_rows[0]) if manifest_rows else ["rank", "root_entry"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(manifest_rows)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "input_root": str(args.input_root.resolve()), "angles": str(args.angles.resolve()),
        "angle_source": args.angle_source, "min_angle_mrad": args.min_angle_mrad,
        "score_min": args.score_min, "coherent_fit_only": args.coherent_fit_only,
        "count": len(manifest_rows), "smoothing": "one ROOT TH2::Smooth(k5a)-equivalent pass per layer",
        "candidates": manifest_rows,
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"count": len(manifest_rows), "output_dir": str(args.output_dir.resolve())}, indent=2))


if __name__ == "__main__":
    main()
