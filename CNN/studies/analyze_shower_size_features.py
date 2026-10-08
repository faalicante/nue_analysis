#!/usr/bin/env python3
"""Compare Poisson-robust shower-size features in signal and hard samples."""

from __future__ import annotations

import argparse
import csv
import glob
import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np
from scipy.ndimage import label
from scipy.stats import mannwhitneyu, poisson

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from training.cnn21d.metrics import slopes_to_angles


POISSON_TAIL_PROBABILITY = 1.0e-4
CROP_SIZE = 20


def poisson_hot_threshold(mu: float, alpha: float = POISSON_TAIL_PROBABILITY) -> int:
    """Smallest integer k for which P(Poisson(mu) >= k) <= alpha."""
    k = max(1, int(poisson.isf(alpha, mu)) + 1)
    while poisson.sf(k - 1, mu) > alpha:
        k += 1
    while k > 1 and poisson.sf(k - 2, mu) <= alpha:
        k -= 1
    return k


def size_features(raw: np.ndarray, mu: float) -> dict[str, float | int]:
    values = np.asarray(raw, dtype=np.float64)
    if values.shape != (57, CROP_SIZE, CROP_SIZE):
        raise ValueError(f"expected [57,20,20], found {values.shape}")
    threshold = poisson_hot_threshold(mu)
    hot = values >= threshold
    residual = values - mu
    layer_net = residual.sum(axis=(1, 2))
    active_layer_threshold = 5.0 * np.sqrt(mu * CROP_SIZE * CROP_SIZE)

    # A 26-neighbour component measures spatial/z coherence.  Q_hot below is
    # the primary size variable because it does not fragment steep showers.
    components, component_count = label(hot, structure=np.ones((3, 3, 3), dtype=bool))
    largest_charge = 0.0
    largest_voxels = 0
    largest_layers = 0
    for component_id in range(1, component_count + 1):
        mask = components == component_id
        charge = float(residual[mask].sum())
        if charge > largest_charge:
            largest_charge = charge
            largest_voxels = int(mask.sum())
            largest_layers = int(np.any(mask, axis=(1, 2)).sum())

    return {
        "background_mu": float(mu),
        "poisson_hot_count_threshold": threshold,
        "net_excess": float(residual.sum()),
        "q_hot": float(residual[hot].sum()),
        "hot_voxels": int(hot.sum()),
        "active_layers_5sigma": int((layer_net > active_layer_threshold).sum()),
        "largest_component_charge": largest_charge,
        "largest_component_voxels": largest_voxels,
        "largest_component_layers": largest_layers,
    }


def background_candidate_rows(volume_glob: str, prediction_glob: str) -> list[dict[str, object]]:
    sources = {path.name: path for path in map(Path, glob.glob(volume_glob))}
    rows: list[dict[str, object]] = []
    for prediction_path in sorted(map(Path, glob.glob(prediction_glob))):
        source_name = prediction_path.name.removesuffix("_t098.h5") + ".h5"
        with h5py.File(sources[source_name], "r") as source, h5py.File(prediction_path, "r") as pred:
            group = pred["candidates"]
            x_positions = np.asarray(pred["grid_x_start_bin"], dtype=int)
            y_positions = np.asarray(pred["grid_y_start_bin"], dtype=int)
            for index in range(len(group["event_id"])):
                event_index = int(group["event_index"][index])
                grid_row = int(group["grid_row"][index])
                grid_col = int(group["grid_col"][index])
                x0 = int(x_positions[grid_col])
                y0 = int(y_positions[grid_row])
                raw = source["volumes_raw"][event_index, :, y0:y0 + 20, x0:x0 + 20]
                row: dict[str, object] = {
                    "sample": "background_scan_candidate",
                    "event_id": int(group["event_id"][index]),
                    "cell_x": int(source["cell_x"][event_index]),
                    "cell_y": int(source["cell_y"][event_index]),
                    "candidate_id": int(group["candidate_id"][index]),
                    "presence_score": float(group["presence_score"][index]),
                    "theta_true_mrad": "",
                    "theta_pred_mrad": float(group["theta_pred_mrad"][index]),
                    "source_entry": event_index,
                }
                row.update(size_features(raw, float(source["background_mu"][event_index])))
                rows.append(row)
    return rows


def signal_rows(volumes: Path, predictions: Path, theta_min: float) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with h5py.File(volumes, "r") as source, h5py.File(predictions, "r") as pred:
        group = pred["candidates"]
        theta_true = slopes_to_angles(np.asarray(source["truth/slope_xy_bins_per_z"]))[0]
        x_positions = np.asarray(pred["grid_x_start_bin"], dtype=int)
        y_positions = np.asarray(pred["grid_y_start_bin"], dtype=int)
        best_by_event: dict[int, int] = {}
        for candidate_index in range(len(group["event_id"])):
            event_index = int(group["event_index"][candidate_index])
            if not bool(group["truth_matched"][candidate_index]) or theta_true[event_index] <= theta_min:
                continue
            previous = best_by_event.get(event_index)
            if previous is None or group["presence_score"][candidate_index] > group["presence_score"][previous]:
                best_by_event[event_index] = candidate_index

        for event_index, candidate_index in best_by_event.items():
            grid_row = int(group["grid_row"][candidate_index])
            grid_col = int(group["grid_col"][candidate_index])
            x0 = int(x_positions[grid_col])
            y0 = int(y_positions[grid_row])
            raw = source["volumes_raw"][event_index, :, y0:y0 + 20, x0:x0 + 20]
            row: dict[str, object] = {
                "sample": "signal_truth_matched",
                "event_id": int(source["event_id"][event_index]),
                "cell_x": "",
                "cell_y": "",
                "candidate_id": int(group["candidate_id"][candidate_index]),
                "presence_score": float(group["presence_score"][candidate_index]),
                "theta_true_mrad": float(theta_true[event_index]),
                "theta_pred_mrad": float(group["theta_pred_mrad"][candidate_index]),
                "source_entry": event_index,
            }
            row.update(size_features(raw, float(source["background_mu"][event_index])))
            rows.append(row)
    return rows


def original_hard_rows(dataset: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with h5py.File(dataset, "r") as source:
        indices = np.flatnonzero(np.asarray(source["sample_type"]) == 1)
        for index in indices:
            raw = source["regions_raw"][index, :, 6:26, 6:26]
            row: dict[str, object] = {
                "sample": "original_hard",
                "event_id": int(source["signal_event_id"][index]),
                "cell_x": "",
                "cell_y": "",
                "candidate_id": int(source["crop_id_in_tile"][index]),
                "presence_score": "",
                "theta_true_mrad": "",
                "theta_pred_mrad": "",
                "source_entry": int(source["source_entry"][index]),
            }
            row.update(size_features(raw, float(source["background_mu"][index, 0])))
            rows.append(row)
    return rows


def numeric(rows: list[dict[str, object]], key: str) -> np.ndarray:
    return np.asarray([float(row[key]) for row in rows], dtype=np.float64)


def summarize_sample(rows: list[dict[str, object]]) -> dict[str, object]:
    summary: dict[str, object] = {"count": len(rows)}
    for feature in (
        "net_excess", "q_hot", "hot_voxels", "active_layers_5sigma",
        "largest_component_charge", "largest_component_voxels", "largest_component_layers",
    ):
        values = numeric(rows, feature)
        summary[feature] = {
            "median": float(np.median(values)),
            "q10": float(np.quantile(values, 0.10)),
            "q25": float(np.quantile(values, 0.25)),
            "q75": float(np.quantile(values, 0.75)),
            "q90": float(np.quantile(values, 0.90)),
        }
    return summary


def auc_probability(signal: np.ndarray, comparison: np.ndarray) -> float:
    statistic = mannwhitneyu(signal, comparison, alternative="greater").statistic
    return float(statistic / (len(signal) * len(comparison)))


def save_plot(samples: dict[str, list[dict[str, object]]], output: Path) -> None:
    colors = {
        "signal": "tab:blue",
        "background_candidates": "tab:red",
        "original_hard": "tab:orange",
    }
    labels = {
        "signal": "signal rilevato, theta MC > 10 mrad",
        "background_candidates": "42 candidati fondo scan",
        "original_hard": "hard originale",
    }
    figure, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    features = (
        ("largest_component_voxels", "voxel hot nella componente maggiore"),
        ("q_hot", "Q_hot [conteggi sopra il Poisson]"),
        ("net_excess", "eccesso netto totale [conteggi]"),
    )
    for axis, (feature, xlabel) in zip(axes, features, strict=True):
        all_values = np.concatenate([numeric(rows, feature) for rows in samples.values()])
        if feature == "net_excess":
            low, high = np.quantile(all_values, [0.001, 0.999])
            bins = np.linspace(low, high, 42)
        else:
            positive = np.maximum(all_values, 0.5)
            low = max(float(np.quantile(positive, 0.001)), 0.5)
            high = float(np.quantile(positive, 0.999))
            bins = np.geomspace(low, high, 42)
        for name, rows in samples.items():
            values = numeric(rows, feature)
            if feature != "net_excess":
                values = np.maximum(values, 0.5)
            axis.hist(
                values, bins=bins, density=True, histtype="step", linewidth=1.8,
                color=colors[name], label=f"{labels[name]} (N={len(rows)})",
            )
        if feature != "net_excess":
            axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel(xlabel)
        axis.set_ylabel("densità normalizzata")
        axis.grid(alpha=0.18)
    axes[0].axvline(8.0, color="black", linestyle="--", linewidth=1.2, label="N componente = 8")
    handles, legend_labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles, legend_labels, loc="upper center", bbox_to_anchor=(0.5, 0.91),
        ncol=2, frameon=False,
    )
    figure.suptitle("Dimensione degli sciami nei crop 20x20 selezionati", fontsize=14, y=0.99)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.80))
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--background-volume-glob",
        default="/Users/fabioali/cernbox/CNN/background_scan_volumes/background_scan_*.h5",
    )
    parser.add_argument(
        "--background-prediction-glob",
        default=("runs/cnn21d_signal_t5_center40/background_scan_full_t098/"
                 "background_scan_*_t098.h5"),
    )
    parser.add_argument("--signal-volumes", type=Path, default=Path("/private/tmp/crop_bw_scan_splits/signal_scan_test.h5"))
    parser.add_argument("--signal-predictions", type=Path, default=Path("runs/cnn21d_signal_t5_center40/scan_volumes_test_gt10_t098.h5"))
    parser.add_argument("--training-dataset", type=Path, default=Path("cnn_dataset_signal_t5_center40.h5"))
    parser.add_argument("--output-dir", type=Path, default=Path("runs/cnn21d_signal_t5_center40/shower_size_comparison"))
    parser.add_argument("--theta-min", type=float, default=10.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    background = background_candidate_rows(args.background_volume_glob, args.background_prediction_glob)
    signal = signal_rows(args.signal_volumes, args.signal_predictions, args.theta_min)
    hard = original_hard_rows(args.training_dataset)
    samples = {"signal": signal, "background_candidates": background, "original_hard": hard}
    all_rows = signal + background + hard

    csv_path = args.output_dir / "shower_size_features.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)

    ranked_path = args.output_dir / "background_candidates_ranked_by_size.csv"
    ranked = sorted(
        background,
        key=lambda row: (float(row["largest_component_voxels"]), float(row["q_hot"])),
        reverse=True,
    )
    with ranked_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["size_rank", *list(ranked[0])])
        writer.writeheader()
        for rank, row in enumerate(ranked, start=1):
            writer.writerow({"size_rank": rank} | row)

    signal_q = numeric(signal, "q_hot")
    background_q = numeric(background, "q_hot")
    hard_q = numeric(hard, "q_hot")
    thresholds = [100, 150, 200, 250, 300, 400, 500]
    threshold_rows = []
    for threshold in thresholds:
        threshold_rows.append({
            "q_hot_threshold": threshold,
            "background_candidates_kept": int((background_q >= threshold).sum()),
            "background_candidates_fraction": float((background_q >= threshold).mean()),
            "signal_fraction": float((signal_q >= threshold).mean()),
            "original_hard_fraction": float((hard_q >= threshold).mean()),
        })

    signal_component = numeric(signal, "largest_component_voxels")
    background_component = numeric(background, "largest_component_voxels")
    hard_component = numeric(hard, "largest_component_voxels")
    component_threshold_rows = []
    for threshold in [3, 5, 8, 10, 12, 15, 20, 25]:
        component_threshold_rows.append({
            "largest_component_voxels_threshold": threshold,
            "background_candidates_kept": int((background_component >= threshold).sum()),
            "background_candidates_fraction": float((background_component >= threshold).mean()),
            "signal_fraction": float((signal_component >= threshold).mean()),
            "original_hard_fraction": float((hard_component >= threshold).mean()),
        })

    report = {
        "feature_definition": {
            "primary": (
                "N_hot_component = number of Poisson-incompatible voxels in the largest "
                "26-connected component"
            ),
            "secondary": "Q_hot = sum(H-mu) over voxels satisfying P(Poisson(mu)>=H) <= 1e-4",
            "poisson_tail_probability": POISSON_TAIL_PROBABILITY,
            "crop_shape": [57, 20, 20],
            "largest_component_connectivity": "26-neighbour in z,y,x among hot voxels",
        },
        "samples": {name: summarize_sample(rows) for name, rows in samples.items()},
        "separation_probability_auc": {
            "q_hot_signal_gt_background_candidates": auc_probability(signal_q, background_q),
            "q_hot_signal_gt_original_hard": auc_probability(signal_q, hard_q),
            "largest_component_voxels_signal_gt_background_candidates": auc_probability(
                numeric(signal, "largest_component_voxels"),
                numeric(background, "largest_component_voxels"),
            ),
        },
        "largest_component_operating_points": component_threshold_rows,
        "q_hot_operating_points": threshold_rows,
        "outputs": {"all_features_csv": str(csv_path), "ranked_background_csv": str(ranked_path)},
    }
    report_path = args.output_dir / "summary.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    plot_path = args.output_dir / "shower_size_distributions.png"
    save_plot(samples, plot_path)
    print(json.dumps({
        "background_candidates": len(background), "signal": len(signal), "original_hard": len(hard),
        "plot": str(plot_path.resolve()), "summary": str(report_path.resolve()),
        "ranked_background": str(ranked_path.resolve()),
    }, indent=2))


if __name__ == "__main__":
    main()
