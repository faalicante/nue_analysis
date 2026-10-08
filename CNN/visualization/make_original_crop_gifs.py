#!/usr/bin/env python3
"""Create smoothed 57-layer GIFs directly from original 32x32 ROOT crops."""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import matplotlib
import numpy as np
import uproot
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from visualization.make_background_candidate_gifs import root_th2_smooth


def render_frames(values: np.ndarray, title: str) -> list[Image.Image]:
    finite = values[np.isfinite(values)]
    vmin = float(np.quantile(finite, 0.02))
    vmax = max(float(np.quantile(finite, 0.998)), vmin + 1.0)
    figure, axis = plt.subplots(figsize=(6.2, 5.9), constrained_layout=True)
    image = axis.imshow(
        values[0], origin="lower", cmap="viridis", vmin=vmin, vmax=vmax,
        extent=(0, 32, 0, 32), interpolation="nearest", aspect="equal",
    )
    colorbar = figure.colorbar(image, ax=axis, shrink=0.84)
    colorbar.set_label("conteggi smussati")
    axis.set_xlabel("x [bin]")
    axis.set_ylabel("y [bin]")
    layer_title = axis.set_title(f"{title}\nlayer 01/57", fontsize=11)
    frames: list[Image.Image] = []
    for z_index in range(values.shape[0]):
        image.set_data(values[z_index])
        layer_title.set_text(f"{title}\nlayer {z_index + 1:02d}/57")
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=100)
        buffer.seek(0)
        frames.append(
            Image.open(buffer).convert("P", palette=Image.Palette.ADAPTIVE, colors=256)
        )
        buffer.close()
    plt.close(figure)
    return frames


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--class-name", choices=("hard", "poisson"), required=True)
    parser.add_argument("--seed", type=int, default=20260831)
    parser.add_argument("--frame-duration-ms", type=int, default=160)
    args = parser.parse_args()
    if args.count <= 0:
        raise ValueError("count must be positive")

    root_file = uproot.open(args.input)
    tree = root_file["samples"]
    entries = int(tree.num_entries)
    if args.count > entries:
        raise ValueError(f"requested {args.count} examples from only {entries}")
    rng = np.random.default_rng(args.seed)
    selected = np.sort(rng.choice(entries, size=args.count, replace=False))

    expected_type = 1 if args.class_name == "hard" else 0
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for ordinal, entry in enumerate(selected, start=1):
        arrays = tree.arrays(
            ["counts", "background_mu", "sample_type", "cell_id", "tag_cell_id"],
            entry_start=int(entry), entry_stop=int(entry) + 1, library="np",
        )
        if int(arrays["sample_type"][0]) != expected_type:
            raise ValueError(f"ROOT entry {entry} does not match requested sample type")
        raw = np.asarray(arrays["counts"][0], dtype=np.float32)
        if raw.shape != (57, 32, 32):
            raise ValueError(f"entry {entry} has shape {raw.shape}, expected (57,32,32)")
        smoothed = root_th2_smooth(raw)
        cell_id = int(arrays["cell_id"][0])
        tag_cell_id = int(arrays["tag_cell_id"][0])
        title = (
            f"{args.class_name} originale 32×32 · ROOT entry {entry}\n"
            f"cell={cell_id} · tag_cell={tag_cell_id} · μ={float(arrays['background_mu'][0]):.3f}"
        )
        frames = render_frames(smoothed, title)
        filename = (
            f"{args.class_name}_{ordinal:02d}_rootentry_{entry:05d}_"
            f"cell_{cell_id:03d}_tag_{tag_cell_id:03d}.gif"
        )
        output = args.output_dir / filename
        durations = [1000] + [args.frame_duration_ms] * 55 + [1000]
        frames[0].save(
            output, save_all=True, append_images=frames[1:], duration=durations,
            loop=0, optimize=False, disposal=2,
        )
        rows.append({
            "ordinal": ordinal,
            "class_name": args.class_name,
            "root_entry": int(entry),
            "sample_type": expected_type,
            "cell_id": cell_id,
            "tag_cell_id": tag_cell_id,
            "background_mu": float(arrays["background_mu"][0]),
            "input_root": str(args.input.resolve()),
            "gif": str(output.resolve()),
        })
        print(json.dumps({"written": ordinal, "root_entry": int(entry), "gif": str(output)}), flush=True)

    root_file.close()

    manifest_csv = args.output_dir / "manifest.csv"
    with manifest_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest_json = args.output_dir / "manifest.json"
    manifest_json.write_text(json.dumps({
        "input_root": str(args.input.resolve()),
        "selection": "uniform random ROOT entries without replacement",
        "seed": args.seed,
        "class_name": args.class_name,
        "count": len(rows),
        "crop_shape": [57, 32, 32],
        "smoothing": "one ROOT TH2::Smooth(k5a)-equivalent pass per layer",
        "tagging_overlays": False,
        "samples": rows,
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "count": len(rows),
        "output_dir": str(args.output_dir.resolve()),
        "manifest_csv": str(manifest_csv.resolve()),
    }, indent=2))


if __name__ == "__main__":
    main()
