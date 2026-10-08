#!/usr/bin/env python3
"""Compare old/new ROOT crops for the previously problematic signal events."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

from studies.fit_signal_centroids import fit_axis_ransac


EVENT_IDS = (5461, 3736, 938, 9255, 5022, 6008)
OLD_ROOT = Path("/Users/fabioali/cernbox/CNN/signal.root")
NEW_ROOT = Path("/Users/fabioali/cernbox/CNN/signal_t5.root")
OUTPUT_DIR = Path("runs/cnn21d_signal_t5_center40/crop_comparison")


def index_by_event(tree: uproot.TTree) -> dict[int, int]:
    event_ids = tree["signal_event_id"].array(library="np")
    return {int(event): int(index) for index, event in enumerate(event_ids)}


def read_entry(tree: uproot.TTree, entry: int) -> dict[str, object]:
    names = ["counts", "background_mu", "slope_x", "slope_y", "zScale"]
    arrays = tree.arrays(names, entry_start=entry, entry_stop=entry + 1, library="np")
    raw = np.asarray(arrays["counts"][0])
    background = np.full(raw.shape[0], float(arrays["background_mu"][0]))
    fit, significance = fit_axis_ransac(raw, background, iterations=10_000, seed=entry)
    return {
        "raw": raw,
        "significance": significance,
        "slope": np.asarray([arrays["slope_x"][0], arrays["slope_y"][0]], dtype=np.float64) * 0.027,
        "z_scale": float(arrays["zScale"][0]),
        "fit": fit,
    }


def draw_projection(axis: plt.Axes, item: dict[str, object], coordinate: str) -> None:
    significance = np.asarray(item["significance"])
    if coordinate == "x":
        projection = significance.max(axis=1).T
        coordinate_size = significance.shape[2]
        component = 0
    else:
        projection = significance.max(axis=2).T
        coordinate_size = significance.shape[1]
        component = 1
    axis.imshow(projection, origin="lower", aspect="auto", cmap="magma", extent=(0, 57, 0, coordinate_size))
    fit = item["fit"]
    z = np.arange(57)
    if fit is not None:
        fit_slope = fit.slope_x if component == 0 else fit.slope_y
        fit_intercept = fit.intercept_x if component == 0 else fit.intercept_y
        axis.plot(z, fit_intercept + fit_slope * z, color="lime", lw=2, label="fit cieco")
        label_slope = np.asarray(item["slope"])[component]
        centre_z = 28.0
        label_intercept = fit_intercept + (fit_slope - label_slope) * centre_z
        axis.plot(z, label_intercept + label_slope * z, "--", color="cyan", lw=2, label="slope label")
    axis.set(xlim=(0, 56), ylim=(0, coordinate_size), xlabel="layer z", ylabel=f"{coordinate} [bin]")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    trees = {"old": uproot.open(OLD_ROOT)["samples"], "new": uproot.open(NEW_ROOT)["samples"]}
    indices = {name: index_by_event(tree) for name, tree in trees.items()}
    results = []
    for event in EVENT_IDS:
        if any(event not in mapping for mapping in indices.values()):
            results.append({"event": event, "missing": [name for name, mapping in indices.items() if event not in mapping]})
            continue
        items = {name: read_entry(trees[name], indices[name][event]) for name in ("old", "new")}
        old_raw = np.asarray(items["old"]["raw"])
        new_raw = np.asarray(items["new"]["raw"])
        changed = not np.array_equal(old_raw, new_raw)

        figure, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True, sharex=True, sharey="col")
        for row, name in enumerate(("old", "new")):
            for column, coordinate in enumerate(("x", "y")):
                draw_projection(axes[row, column], items[name], coordinate)
            fit = items[name]["fit"]
            fit_text = "fit assente" if fit is None else (
                f"fit=({fit.slope_x:+.3f},{fit.slope_y:+.3f}), θ={fit.theta_mrad:.1f} mrad, "
                f"slice={fit.slices}, reliable={fit.reliable}"
            )
            axes[row, 0].set_title(f"{name}: z–x | {fit_text}")
            axes[row, 1].set_title(f"{name}: z–y | zScale={items[name]['z_scale']:.2f}")
        axes[0, 0].legend(fontsize=8)
        figure.suptitle(f"Signal event {event} — crop cambiato: {changed}")
        path = OUTPUT_DIR / f"event_{event}_old_vs_center40.png"
        figure.savefig(path, dpi=170)
        plt.close(figure)

        event_result: dict[str, object] = {
            "event": event,
            "old_entry": indices["old"][event],
            "new_entry": indices["new"][event],
            "counts_identical": not changed,
            "changed_voxels": int(np.count_nonzero(old_raw != new_raw)),
            "l1_count_difference": int(np.abs(old_raw.astype(np.int64) - new_raw.astype(np.int64)).sum()),
            "plot": str(path),
        }
        for name in ("old", "new"):
            fit = items[name]["fit"]
            event_result[name] = {
                "z_scale": items[name]["z_scale"],
                "fit": None if fit is None else asdict(fit),
            }
        results.append(event_result)
    output = OUTPUT_DIR / "comparison.json"
    output.write_text(json.dumps({"events": results}, indent=2) + "\n")
    print(json.dumps({"output": str(output), "events": results}, indent=2))


if __name__ == "__main__":
    main()
