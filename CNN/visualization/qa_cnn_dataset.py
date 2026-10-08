#!/usr/bin/env python3
"""Visual QA for central versus training-augmented 2+1D CNN inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np

from training.cnn_dataset import CNNDataConfig, HDF5CNN21DDataset


CLASS_NAMES = {0: "poisson", 1: "hard", 2: "signal"}


def _choose_stratified_indices(
    dataset: HDF5CNN21DDataset, count: int, seed: int
) -> list[int]:
    rng = np.random.default_rng(seed)
    with h5py.File(dataset.path, "r") as hdf5:
        sample_types = np.asarray(hdf5["sample_type"][dataset.indices])
    pools = {
        value: rng.permutation(np.flatnonzero(sample_types == value)).tolist()
        for value in CLASS_NAMES
    }
    selected: list[int] = []
    while len(selected) < min(count, len(dataset)):
        made_progress = False
        for value in CLASS_NAMES:
            if pools[value] and len(selected) < count:
                selected.append(int(pools[value].pop()))
                made_progress = True
        if not made_progress:
            break
    return selected


def _normalization_for(before: np.ndarray, after: np.ndarray) -> Normalize:
    combined = np.concatenate((before.ravel(), after.ravel()))
    low, high = np.quantile(combined, [0.01, 0.99])
    if not np.isfinite(low) or not np.isfinite(high) or low == high:
        low, high = float(np.min(combined)), float(np.max(combined) + 1.0)
    if low < 0 < high:
        extent = max(abs(float(low)), abs(float(high)))
        return TwoSlopeNorm(vmin=-extent, vcenter=0.0, vmax=extent)
    return Normalize(vmin=float(low), vmax=float(high))


def _draw_slope(axis: plt.Axes, slope_xy: np.ndarray) -> None:
    """Draw displacement over 20 z bins, capped only to stay inside the panel."""

    displacement = np.asarray(slope_xy, dtype=np.float64) * 20.0
    magnitude = float(np.linalg.norm(displacement))
    if magnitude == 0:
        axis.text(
            0.03,
            0.04,
            "slope = (0, 0)",
            transform=axis.transAxes,
            color="white",
            fontsize=9,
            bbox={"facecolor": "black", "alpha": 0.55, "pad": 2},
        )
        return
    capped = displacement * min(1.0, 7.0 / magnitude)
    axis.arrow(
        9.5,
        9.5,
        capped[0],
        capped[1],
        width=0.12,
        head_width=0.75,
        head_length=0.8,
        length_includes_head=True,
        facecolor="cyan",
        edgecolor="black",
        linewidth=0.8,
        zorder=4,
    )
    cap_note = ", display capped" if magnitude > 7.0 else ""
    axis.text(
        0.03,
        0.04,
        f"slope=({slope_xy[0]:+.3f}, {slope_xy[1]:+.3f})\n"
        f"arrow = 20 z bins{cap_note}",
        transform=axis.transAxes,
        color="white",
        fontsize=8,
        bbox={"facecolor": "black", "alpha": 0.55, "pad": 2},
    )


def render_example(
    before: dict[str, Any], after: dict[str, Any], output: Path
) -> None:
    before_xy = before["volume"][0].sum(dim=0).numpy()
    after_xy = after["volume"][0].sum(dim=0).numpy()
    norm = _normalization_for(before_xy, after_xy)
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 5.2), constrained_layout=True)
    images = []
    for axis, projection, title, slope in (
        (
            axes[0],
            before_xy,
            "Prima: crop centrale, nessuna augmentation",
            before["slope_xy"].numpy(),
        ),
        (
            axes[1],
            after_xy,
            "Dopo: crop traslato + trasformazione xy",
            after["slope_xy"].numpy(),
        ),
    ):
        image = axis.imshow(
            projection,
            origin="lower",
            interpolation="nearest",
            cmap="coolwarm",
            norm=norm,
            extent=(-0.5, 19.5, -0.5, 19.5),
        )
        images.append(image)
        axis.set_xlabel("x bin (ROOT)")
        axis.set_ylabel("y bin (ROOT)")
        axis.set_title(title, fontsize=10)
        _draw_slope(axis, slope)

    metadata = after["metadata"]
    sample_type = int(metadata["sample_type"])
    transform = (
        f"offset xy={metadata['crop_offset_xy'].tolist()} | "
        f"rot={int(metadata['rotation_quarter_turns_ccw']) * 90}° CCW | "
        f"reflect x/y={bool(metadata['reflected_x'])}/{bool(metadata['reflected_y'])}"
    )
    visibility = (
        f"visible fraction={float(metadata['positive_visible_fraction']):.3f}, "
        f"warning={bool(metadata['positive_may_be_insufficiently_visible'])}"
    )
    figure.suptitle(
        f"{CLASS_NAMES.get(sample_type, str(sample_type))} | "
        f"HDF5 row {int(metadata['hdf5_index'])}\n{transform}\n{visibility}",
        fontsize=9,
    )
    figure.colorbar(images[-1], ax=axes, shrink=0.88, label="Σz normalized slice values")
    figure.savefig(output, dpi=160)
    plt.close(figure)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hdf5", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("cnn_dataloader_qa"))
    parser.add_argument("--count", type=int, default=6)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--epoch", type=int, default=0)
    parser.add_argument("--max-jitter", type=int, default=5)
    parser.add_argument(
        "--positive-max-jitter",
        type=int,
        default=None,
        help="optional reduced jitter for sample_type 1 or 2 (for example 3)",
    )
    parser.add_argument("--min-visible-fraction", type=float, default=0.70)
    parser.add_argument(
        "--skip-visibility-audit",
        action="store_true",
        help="skip the split-wide worst-case positive visibility count",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.count <= 0:
        raise SystemExit("--count must be positive")
    config = CNNDataConfig(
        max_jitter=args.max_jitter,
        positive_max_jitter=args.positive_max_jitter,
        min_positive_visible_fraction=args.min_visible_fraction,
        seed=args.seed,
    )
    central = HDF5CNN21DDataset(
        args.hdf5, split="train", training=False, config=config
    )
    augmented = HDF5CNN21DDataset(
        args.hdf5, split="train", training=True, config=config
    )
    augmented.set_epoch(args.epoch)
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = []
    for ordinal, local_index in enumerate(
        _choose_stratified_indices(augmented, args.count, args.seed)
    ):
        before = central[local_index]
        after = augmented[local_index]
        sample_type = int(after["metadata"]["sample_type"])
        hdf5_index = int(after["metadata"]["hdf5_index"])
        output = output_dir / (
            f"qa_{ordinal:02d}_{CLASS_NAMES.get(sample_type, sample_type)}_"
            f"hdf5_{hdf5_index:06d}.png"
        )
        render_example(before, after, output)
        outputs.append(str(output))

    report = None
    if not args.skip_visibility_audit:
        report = augmented.audit_positive_visibility()
        report_path = output_dir / "visibility_report.json"
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    central.close()
    augmented.close()
    print(json.dumps({"images": outputs, "visibility_report": report}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
