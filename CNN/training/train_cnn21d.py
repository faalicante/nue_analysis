#!/usr/bin/env python3
"""Train and diagnose the two-headed 2+1D CNN.

Modes are intentionally explicit: ``pilot`` is subset-only, while ``full`` is
the later production command. No mode silently changes a configured subset.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import time
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/cnn21d-matplotlib")

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import yaml
from torch import nn
from torch.utils.data import DataLoader, Subset

from training.cnn_dataset import CNNDataConfig, HDF5CNN21DDataset, seed_dataloader_worker
from training.cnn21d.losses import MultiTaskLoss
from training.cnn21d.metrics import (
    NEGATIVE_TYPE_NAMES,
    SAMPLE_TYPE_NAMES,
    classification_metrics,
    curve_points,
    regression_metrics,
)
from training.cnn21d.model import CNN21DConfig, TwoHeadCNN21D


def set_seed(seed: int, deterministic: bool = True) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if deterministic:
        torch.use_deterministic_algorithms(True, warn_only=True)
        torch.backends.cudnn.benchmark = False


def choose_device(requested: str) -> torch.device:
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError("YAML root must be a mapping")
    return config


def build_model(config: dict[str, Any]) -> TwoHeadCNN21D:
    model = config["model"]
    parsed = CNN21DConfig(
        in_channels=int(model.get("in_channels", 1)),
        channels=tuple(int(x) for x in model["channels"]),
        group_norm_groups=tuple(int(x) for x in model["group_norm_groups"]),
        pool_kernels=tuple(tuple(int(y) for y in x) for x in model["pool_kernels"]),
        head_hidden=int(model.get("head_hidden", 64)),
        dropout=float(model.get("dropout", 0.0)),
    )
    return TwoHeadCNN21D(parsed)


def _stratified_local_indices(
    dataset: HDF5CNN21DDataset,
    size: int | None,
    seed: int,
    include_sample_types: tuple[int, ...],
    positive_sample_types: tuple[int, ...],
    signal_theta_range_mrad: tuple[float, float] | None,
    sample_type_fractions: dict[int, float] | None = None,
) -> list[int]:
    with h5py.File(dataset.path, "r") as hdf5:
        sample_types = np.asarray(hdf5["sample_type"][dataset.indices], dtype=np.int64)
        slope_xy = np.asarray(hdf5["slope_xy"][dataset.indices], dtype=np.float64)
    eligible_mask = np.isin(sample_types, include_sample_types)
    if signal_theta_range_mrad is not None:
        minimum, maximum = signal_theta_range_mrad
        theta_mrad = 1000.0 * np.arctan(np.linalg.norm(slope_xy, axis=1) / 27.0)
        positive = np.isin(sample_types, positive_sample_types)
        eligible_mask &= ~positive | ((theta_mrad > minimum) & (theta_mrad < maximum))
    eligible = np.flatnonzero(eligible_mask)
    if eligible.size == 0:
        raise ValueError(f"no samples remain after selecting sample_type in {include_sample_types}")
    if size is None or size >= eligible.size:
        return eligible.tolist()
    labels = np.isin(sample_types, positive_sample_types).astype(np.int64)
    # Keeping sample_type in the stratum permits an explicit small Poisson veto
    # quota without allowing it to dominate the physically relevant hard class.
    strata = labels * 10 + sample_types
    rng = np.random.default_rng(seed)
    groups = {int(code): eligible[strata[eligible] == code] for code in np.unique(strata[eligible])}
    if sample_type_fractions is None:
        allocations = {
            code: max(1, int(round(size * len(values) / eligible.size)))
            for code, values in groups.items()
        }
    else:
        if any(value < 0 for value in sample_type_fractions.values()):
            raise ValueError("sample_type_fractions must be non-negative")
        total_fraction = sum(sample_type_fractions.values())
        if not np.isclose(total_fraction, 1.0):
            raise ValueError("sample_type_fractions must sum to one")
        allocations = {}
        for code, values in groups.items():
            sample_type_code = int(code % 10)
            fraction = sample_type_fractions.get(sample_type_code, 0.0)
            allocations[code] = min(len(values), int(round(size * fraction)))
        if not any(allocations.values()):
            raise ValueError("sample_type_fractions allocate no eligible examples")
    while sum(allocations.values()) > size:
        code = max((c for c in allocations if allocations[c] > 1), key=lambda c: allocations[c])
        allocations[code] -= 1
    while sum(allocations.values()) < size:
        code = max(groups, key=lambda c: len(groups[c]) - allocations[c])
        allocations[code] += 1
    selected: list[int] = []
    for code, values in groups.items():
        selected.extend(rng.choice(values, size=min(allocations[code], len(values)), replace=False).tolist())
    rng.shuffle(selected)
    return selected


def make_loader(
    path: Path,
    split: str,
    data_cfg: CNNDataConfig,
    *,
    batch_size: int,
    workers: int,
    seed: int,
    subset_size: int | None,
    training: bool,
    shuffle: bool,
    include_sample_types: tuple[int, ...],
    positive_sample_types: tuple[int, ...],
    signal_theta_range_mrad: tuple[float, float] | None,
    sample_type_fractions: dict[int, float] | None = None,
) -> DataLoader:
    dataset = HDF5CNN21DDataset(path, split=split, training=training, config=data_cfg)  # type: ignore[arg-type]
    indices = _stratified_local_indices(
        dataset, subset_size, seed, include_sample_types, positive_sample_types,
        signal_theta_range_mrad, sample_type_fractions,
    )
    generator = torch.Generator().manual_seed(seed)
    return DataLoader(
        Subset(dataset, indices),
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=workers,
        pin_memory=False,
        worker_init_fn=seed_dataloader_worker,
        generator=generator,
    )


def sample_types_to_labels(sample_type: torch.Tensor, positive_sample_types: tuple[int, ...]) -> torch.Tensor:
    label = torch.zeros_like(sample_type, dtype=torch.bool)
    for value in positive_sample_types:
        label |= sample_type == value
    return label.to(dtype=torch.float32)


def build_task_targets(
    sample_type: torch.Tensor,
    slope_xy: torch.Tensor,
    positive_sample_types: tuple[int, ...],
    presence_positive_theta_min_mrad: float | None = None,
    presence_low_angle_policy: str = "ignore",
    presence_negative_theta_max_mrad: float | None = None,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Build BCE labels/mask and the independent signal regression mask.

    When a theta minimum is configured, signal examples below it always remain
    valid regression examples. For BCE they can either be ignored (legacy
    behaviour) or explicitly labelled negative. Setting a negative-theta
    maximum together with the negative policy creates an ignored transition
    band between that maximum and the positive-theta minimum.
    """
    if presence_low_angle_policy not in {"ignore", "negative"}:
        raise ValueError("presence_low_angle_policy must be 'ignore' or 'negative'")
    if presence_negative_theta_max_mrad is not None:
        if presence_positive_theta_min_mrad is None:
            raise ValueError(
                "presence_negative_theta_max_mrad requires "
                "presence_positive_theta_min_mrad"
            )
        if not 0.0 <= presence_negative_theta_max_mrad <= presence_positive_theta_min_mrad:
            raise ValueError(
                "presence_negative_theta_max_mrad must lie between zero and "
                "presence_positive_theta_min_mrad"
            )
        if presence_low_angle_policy != "negative":
            raise ValueError(
                "presence_negative_theta_max_mrad requires "
                "presence_low_angle_policy='negative'"
            )
    signal_mask = sample_types_to_labels(sample_type, positive_sample_types).bool()
    regression_mask = signal_mask.clone()
    if presence_positive_theta_min_mrad is None:
        return signal_mask.float(), torch.ones_like(signal_mask), regression_mask
    theta_mrad = 1000.0 * torch.atan(torch.linalg.vector_norm(slope_xy, dim=-1) / 27.0)
    high_angle = signal_mask & (theta_mrad > float(presence_positive_theta_min_mrad))
    hard_mask = ~signal_mask
    if presence_low_angle_policy == "ignore":
        presence_mask = hard_mask | high_angle
    elif presence_negative_theta_max_mrad is None:
        presence_mask = torch.ones_like(signal_mask)
    else:
        low_angle = signal_mask & (theta_mrad < presence_negative_theta_max_mrad)
        presence_mask = hard_mask | low_angle | high_angle
    return high_angle.float(), presence_mask, regression_mask


def move_batch(
    batch: dict[str, Any], device: torch.device, positive_sample_types: tuple[int, ...],
    presence_positive_theta_min_mrad: float | None = None,
    presence_low_angle_policy: str = "ignore",
    presence_negative_theta_max_mrad: float | None = None,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    sample_type = batch.get("sample_type", batch["metadata"]["sample_type"])
    slope = batch["slope_xy"]
    labels, presence_mask, regression_mask = build_task_targets(
        sample_type, slope, positive_sample_types, presence_positive_theta_min_mrad,
        presence_low_angle_policy, presence_negative_theta_max_mrad,
    )
    return (
        batch["volume"].to(device),
        labels.to(device),
        slope.to(device),
        presence_mask.to(device),
        regression_mask.to(device),
    )


def train_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: MultiTaskLoss,
    device: torch.device,
    epoch: int,
    positive_sample_types: tuple[int, ...],
    presence_positive_theta_min_mrad: float | None = None,
    presence_low_angle_policy: str = "ignore",
    presence_negative_theta_max_mrad: float | None = None,
) -> dict[str, float]:
    model.train()
    base = loader.dataset.dataset if isinstance(loader.dataset, Subset) else loader.dataset
    if hasattr(base, "set_epoch"):
        base.set_epoch(epoch)
    sums = {"total": 0.0, "presence": 0.0, "slope": 0.0}
    count = 0
    for batch in loader:
        volume, presence, slope, presence_mask, regression_mask = move_batch(
            batch, device, positive_sample_types, presence_positive_theta_min_mrad,
            presence_low_angle_policy, presence_negative_theta_max_mrad,
        )
        optimizer.zero_grad(set_to_none=True)
        losses = criterion(
            model(volume), presence, slope,
            presence_mask=presence_mask, regression_mask=regression_mask,
        )
        losses["total"].backward()
        optimizer.step()
        batch_size = volume.shape[0]
        for key in sums:
            sums[key] += float(losses[key].detach().cpu()) * batch_size
        count += batch_size
    return {key: value / max(count, 1) for key, value in sums.items()}


@torch.no_grad()
def predict(
    model: nn.Module,
    loader: DataLoader,
    criterion: MultiTaskLoss,
    device: torch.device,
    *,
    z_transform: str = "original",
    z_seed: int = 12345,
    positive_sample_types: tuple[int, ...],
    presence_positive_theta_min_mrad: float | None = None,
    presence_low_angle_policy: str = "ignore",
    presence_negative_theta_max_mrad: float | None = None,
) -> dict[str, Any]:
    model.eval()
    collected: dict[str, list[np.ndarray]] = {
        key: [] for key in ("labels", "classification_mask", "regression_mask", "probabilities", "slope_truth", "slope_prediction", "negative_type", "sample_type", "border_fraction")
    }
    losses: list[float] = []
    rng = np.random.default_rng(z_seed)
    fixed_permutation: np.ndarray | None = None
    for batch in loader:
        volume, presence, slope, presence_mask, regression_mask = move_batch(
            batch, device, positive_sample_types, presence_positive_theta_min_mrad,
            presence_low_angle_policy, presence_negative_theta_max_mrad,
        )
        if z_transform == "reverse":
            volume = volume.flip(2)
        elif z_transform == "permute":
            if fixed_permutation is None:
                fixed_permutation = rng.permutation(volume.shape[2])
            volume = volume[:, :, torch.as_tensor(fixed_permutation, device=device)]
        elif z_transform in {
            "drop_one_random", "drop_two_isolated", "drop_two_consecutive", "drop_peak"
        }:
            background_mu = batch["metadata"]["background_mu"].to(device)
            missing_value = -torch.sqrt(background_mu).view(-1, 1, 1, 1)
            z_count = volume.shape[2]
            for batch_index in range(volume.shape[0]):
                if z_transform == "drop_one_random":
                    dropped = [int(rng.integers(0, z_count))]
                elif z_transform == "drop_two_isolated":
                    first = int(rng.integers(0, z_count))
                    candidates = np.asarray([
                        z for z in range(z_count) if abs(z - first) > 1
                    ])
                    dropped = [first, int(rng.choice(candidates))]
                elif z_transform == "drop_two_consecutive":
                    first = int(rng.integers(0, z_count - 1))
                    dropped = [first, first + 1]
                else:
                    activation_by_z = volume[batch_index].clamp_min(0).sum(dim=(0, 2, 3))
                    dropped = [int(torch.argmax(activation_by_z).item())]
                volume[batch_index, :, dropped] = missing_value[batch_index]
        elif z_transform != "original":
            raise ValueError(f"unknown z transform: {z_transform}")
        outputs = model(volume)
        losses.append(float(criterion(
            outputs, presence, slope,
            presence_mask=presence_mask, regression_mask=regression_mask,
        )["total"].cpu()))
        probabilities = torch.sigmoid(outputs["presence_logit"])
        # A reproducible edge/central proxy: fraction of positive activation in a
        # 3-pixel xy rim. It does not use labels or model predictions.
        activation = volume.clamp_min(0).sum(dim=(1, 2))
        rim = activation.clone()
        rim[:, 3:-3, 3:-3] = 0
        border_fraction = rim.sum((1, 2)) / activation.sum((1, 2)).clamp_min(1e-12)
        metadata = batch["metadata"]
        arrays = {
            "labels": presence,
            "classification_mask": presence_mask,
            "regression_mask": regression_mask,
            "probabilities": probabilities,
            "slope_truth": slope,
            "slope_prediction": outputs["slope_xy"],
            "negative_type": metadata["negative_type"],
            "sample_type": metadata["sample_type"],
            "border_fraction": border_fraction,
        }
        for key, value in arrays.items():
            collected[key].append(value.detach().cpu().numpy())
    result = {key: np.concatenate(values) for key, values in collected.items()}
    result["loss"] = float(np.mean(losses)) if losses else float("nan")
    return result


def result_metrics(
    result: dict[str, Any], threshold: float, selection: np.ndarray | None = None
) -> dict[str, Any]:
    classification_mask = result.get("classification_mask", np.ones_like(result["labels"], dtype=bool)).astype(bool)
    regression_mask = result.get("regression_mask", result["labels"] > 0.5).astype(bool)
    if selection is not None:
        classification_mask &= selection
        regression_mask &= selection
    return {
        **classification_metrics(
            result["labels"][classification_mask], result["probabilities"][classification_mask], threshold
        ),
        **regression_metrics(
            regression_mask.astype(np.float32), result["slope_truth"], result["slope_prediction"]
        ),
    }


def balanced_theta_mae_mrad(
    result: dict[str, Any], bin_edges_mrad: tuple[float, ...]
) -> tuple[float, list[dict[str, Any]]]:
    """Equal-bin angular MAE, preventing abundant mid-angle events dominating."""

    mask = result.get("regression_mask", result["labels"] > 0.5).astype(bool)
    truth = 1000.0 * np.arctan(np.linalg.norm(result["slope_truth"], axis=1) / 27.0)
    prediction = 1000.0 * np.arctan(
        np.linalg.norm(result["slope_prediction"], axis=1) / 27.0
    )
    rows: list[dict[str, Any]] = []
    values: list[float] = []
    for low, high in zip(bin_edges_mrad, bin_edges_mrad[1:]):
        selected = mask & (truth >= low) & (truth < high)
        mae = float(np.abs(prediction[selected] - truth[selected]).mean()) if selected.any() else None
        rows.append({"low_mrad": low, "high_mrad": high, "count": int(selected.sum()), "mae_mrad": mae})
        if mae is not None:
            values.append(mae)
    return (float(np.mean(values)) if values else float("inf")), rows


def combine_head_predictions(
    classification_result: dict[str, Any], regression_result: dict[str, Any]
) -> dict[str, Any]:
    """Use the AUPRC checkpoint for presence and regression checkpoint for slope."""

    combined = dict(classification_result)
    combined["slope_prediction"] = regression_result["slope_prediction"]
    return combined


def resolve_regression_balance(
    loader: DataLoader,
    positive_sample_types: tuple[int, ...],
    balance_config: dict[str, Any] | None,
) -> dict[str, Any]:
    if not balance_config or not bool(balance_config.get("enabled", False)):
        return {"enabled": False}
    edges = np.asarray(balance_config["theta_bin_edges_mrad"], dtype=np.float64)
    if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("regression balance theta_bin_edges_mrad must increase")
    subset = loader.dataset
    base = subset.dataset if isinstance(subset, Subset) else subset
    local_indices = np.asarray(subset.indices if isinstance(subset, Subset) else np.arange(len(base)))
    raw_indices = base.indices[local_indices]
    # h5py fancy indexing requires monotonically increasing coordinates; the
    # order is irrelevant for histogram counts.
    raw_indices = np.sort(raw_indices)
    with h5py.File(base.path, "r") as hdf5:
        sample_type = np.asarray(hdf5["sample_type"][raw_indices])
        slope = np.asarray(hdf5["slope_xy"][raw_indices], dtype=np.float64)
    signal = np.isin(sample_type, positive_sample_types)
    theta = 1000.0 * np.arctan(np.linalg.norm(slope[signal], axis=1) / 27.0)
    counts, _ = np.histogram(theta, bins=edges)
    raw_weights = np.divide(
        1.0, np.sqrt(counts), out=np.zeros_like(counts, dtype=np.float64), where=counts > 0
    )
    if not np.any(raw_weights):
        raise ValueError("regression balance found no signal examples")
    normalizer = np.sum(counts * raw_weights) / max(np.sum(counts), 1)
    raw_weights /= normalizer
    weights = np.clip(
        raw_weights,
        float(balance_config.get("minimum_weight", 0.5)),
        float(balance_config.get("maximum_weight", 3.0)),
    )
    weights /= np.sum(counts * weights) / max(np.sum(counts), 1)
    return {
        "enabled": True,
        "theta_bin_edges_mrad": edges.tolist(),
        "counts": counts.tolist(),
        "weights": weights.tolist(),
        "scheme": "clipped_inverse_sqrt_frequency",
    }


def grouped_diagnostics(result: dict[str, Any], threshold: float) -> dict[str, Any]:
    labels = result["labels"]
    probabilities = result["probabilities"]
    groups: dict[str, Any] = {"negative_type": {}, "sample_type": {}, "spatial": {}}
    for code in np.unique(result["negative_type"]):
        mask = (result["negative_type"] == code) & result["classification_mask"].astype(bool)
        groups["negative_type"][NEGATIVE_TYPE_NAMES.get(int(code), str(code))] = {
            "count": int(mask.sum()),
            **classification_metrics(labels[mask], probabilities[mask], threshold),
        }
    for code in np.unique(result["sample_type"]):
        mask = result["sample_type"] == code
        groups["sample_type"][SAMPLE_TYPE_NAMES.get(int(code), str(code))] = {
            "count": int(mask.sum()),
            "mean_positive_probability": float(probabilities[mask].mean()),
        }
    positive = result.get("regression_mask", labels > 0.5).astype(bool)
    if positive.any():
        low, high = np.quantile(result["border_fraction"][positive], [0.25, 0.75])
        for name, mask in {
            "central_lowest_border_quartile": positive & (result["border_fraction"] <= low),
            "edge_highest_border_quartile": positive & (result["border_fraction"] >= high),
        }.items():
            groups["spatial"][name] = {
                "count": int(mask.sum()),
                "border_fraction_mean": float(result["border_fraction"][mask].mean()),
                "metrics": result_metrics(result, threshold, selection=mask),
            }
    return groups


def save_plots(output_dir: Path, history: list[dict[str, Any]], result: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if history:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        epochs = [row["epoch"] for row in history]
        for split in ("train", "validation"):
            axes[0].plot(epochs, [row[split]["total"] for row in history], label=split)
        axes[0].set(xlabel="epoch", ylabel="loss", title="Training history")
        axes[0].legend()
        axes[1].plot(epochs, [row["validation_metrics"]["auroc"] for row in history], label="AUROC")
        axes[1].plot(epochs, [row["validation_metrics"]["auprc"] for row in history], label="AUPRC")
        axes[1].set(xlabel="epoch", ylabel="score", ylim=(0, 1), title="Validation")
        axes[1].legend()
        fig.tight_layout(); fig.savefig(output_dir / "training_history.png", dpi=160); plt.close(fig)
    classification_mask = result.get("classification_mask", np.ones_like(result["labels"], dtype=bool)).astype(bool)
    curves = curve_points(result["labels"][classification_mask], result["probabilities"][classification_mask])
    if curves:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        axes[0].plot(curves["fpr"], curves["tpr"]); axes[0].plot([0, 1], [0, 1], "--", color="gray")
        axes[0].set(xlabel="false positive rate", ylabel="true positive rate", title="ROC")
        axes[1].plot(curves["pr_recall"], curves["pr_precision"])
        axes[1].set(xlabel="recall", ylabel="precision", title="Precision-recall")
        fig.tight_layout(); fig.savefig(output_dir / "roc_pr_curves.png", dpi=160); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4))
    for code, name in SAMPLE_TYPE_NAMES.items():
        values = result["probabilities"][result["sample_type"] == code]
        if values.size:
            ax.hist(values, bins=20, range=(0, 1), alpha=0.45, density=True, label=name)
    ax.set(
        xlabel="predicted positive probability P(sample_type=2)",
        ylabel="density",
        title="Hard-background/signal prediction distributions",
    )
    ax.legend(); fig.tight_layout(); fig.savefig(output_dir / "prediction_distributions.png", dpi=160); plt.close(fig)


def run_training(config: dict[str, Any], mode: str) -> dict[str, Any]:
    seed = int(config["seed"])
    set_seed(seed, bool(config.get("deterministic", True)))
    device = choose_device(str(config.get("device", "auto")))
    path = Path(config["data"]["hdf5_path"]).expanduser().resolve()
    output_dir = Path(config["output_dir"]).expanduser().resolve() / mode
    output_dir.mkdir(parents=True, exist_ok=True)
    data_cfg = CNNDataConfig(**config["data"].get("dataset", {}), seed=seed)
    include_sample_types = tuple(int(value) for value in config["task"]["include_sample_types"])
    positive_sample_types = tuple(int(value) for value in config["task"]["positive_sample_types"])
    theta_values = config["task"].get("signal_theta_range_mrad")
    signal_theta_range_mrad = None if theta_values is None else (float(theta_values[0]), float(theta_values[1]))
    presence_theta_value = config["task"].get("presence_positive_theta_min_mrad")
    presence_positive_theta_min_mrad = None if presence_theta_value is None else float(presence_theta_value)
    presence_low_angle_policy = str(
        config["task"].get("presence_low_angle_policy", "ignore")
    )
    negative_theta_value = config["task"].get("presence_negative_theta_max_mrad")
    presence_negative_theta_max_mrad = (
        None if negative_theta_value is None else float(negative_theta_value)
    )
    if presence_low_angle_policy not in {"ignore", "negative"}:
        raise ValueError("task.presence_low_angle_policy must be 'ignore' or 'negative'")
    if presence_negative_theta_max_mrad is not None:
        if presence_positive_theta_min_mrad is None:
            raise ValueError(
                "task.presence_negative_theta_max_mrad requires "
                "task.presence_positive_theta_min_mrad"
            )
        if not 0.0 <= presence_negative_theta_max_mrad <= presence_positive_theta_min_mrad:
            raise ValueError(
                "task.presence_negative_theta_max_mrad must lie between zero and "
                "task.presence_positive_theta_min_mrad"
            )
        if presence_low_angle_policy != "negative":
            raise ValueError(
                "task.presence_negative_theta_max_mrad requires "
                "task.presence_low_angle_policy='negative'"
            )
    if not set(positive_sample_types).issubset(include_sample_types):
        raise ValueError("positive_sample_types must be a subset of include_sample_types")
    train_cfg = config["training"]
    subset_cfg = config["subsets"][mode]
    batch_size = int(train_cfg["batch_size"])
    workers = int(train_cfg.get("num_workers", 0))
    fraction_values = train_cfg.get("train_sample_type_fractions")
    train_sample_type_fractions = None if fraction_values is None else {
        int(key): float(value) for key, value in fraction_values.items()
    }
    loaders = {
        split: make_loader(
            path, split, data_cfg,
            batch_size=batch_size, workers=workers, seed=seed + i,
            subset_size=subset_cfg.get(split), training=(split == "train"), shuffle=(split == "train"),
            include_sample_types=include_sample_types,
            positive_sample_types=positive_sample_types,
            signal_theta_range_mrad=signal_theta_range_mrad,
            sample_type_fractions=(
                train_sample_type_fractions if split == "train" else None
            ),
        )
        for i, split in enumerate(("train", "validation", "test"))
    }
    model = build_model(config).to(device)
    regression_balance = resolve_regression_balance(
        loaders["train"], positive_sample_types, train_cfg.get("regression_balance")
    )
    loss_config = dict(config["loss"])
    if regression_balance["enabled"]:
        loss_config.update({
            "regression_theta_bin_edges_mrad": regression_balance["theta_bin_edges_mrad"],
            "regression_theta_bin_weights": regression_balance["weights"],
        })
    criterion = MultiTaskLoss(**loss_config).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=float(train_cfg["learning_rate"]), weight_decay=float(train_cfg["weight_decay"]))
    threshold = float(config["evaluation"]["threshold"])
    epochs = int(subset_cfg["epochs"])
    patience = int(train_cfg["early_stopping_patience"])
    classification_monitor = str(train_cfg.get("classification_checkpoint_monitor", "auprc"))
    regression_monitor = str(train_cfg.get(
        "regression_checkpoint_monitor", "balanced_theta_mae_mrad"
    ))
    if classification_monitor != "auprc":
        raise ValueError("classification_checkpoint_monitor currently supports only 'auprc'")
    if regression_monitor not in {"balanced_theta_mae_mrad", "theta_mae_mrad"}:
        raise ValueError("unsupported regression_checkpoint_monitor")
    regression_edges = tuple(float(value) for value in (
        regression_balance.get("theta_bin_edges_mrad")
        or train_cfg.get("regression_checkpoint_theta_bin_edges_mrad", [0, 5, 10, 20, 50, 100, 200])
    ))
    # Optional warm-start from a previous interrupted run.  Checkpoints in
    # this project intentionally contain model weights/monitor state only;
    # the AdamW moments are reinitialised, so this is a robust continuation
    # rather than an exact optimizer-state resume.
    resume_path = train_cfg.get("resume_checkpoint")
    resume_epoch = 0
    best_classification_score = -float("inf")
    best_regression_score = float("inf")
    if resume_path:
        resume_ckpt = torch.load(Path(resume_path), map_location=device, weights_only=False)
        model.load_state_dict(resume_ckpt["model_state"])
        resume_epoch = int(resume_ckpt.get("epoch", 0))
        if resume_ckpt.get("checkpoint_monitor") == "auprc":
            best_classification_score = float(resume_ckpt.get("checkpoint_score", best_classification_score))
        print(json.dumps({"resume_checkpoint": str(resume_path), "resume_epoch": resume_epoch}), flush=True)
    min_delta = float(train_cfg.get("early_stopping_min_delta", 0.0))
    classification_stale = 0
    regression_stale = 0
    history: list[dict[str, Any]] = []
    checkpoint_path = output_dir / "best_model.pt"
    classification_checkpoint_path = output_dir / "best_classification_model.pt"
    regression_checkpoint_path = output_dir / "best_regression_model.pt"
    started = time.time()
    for epoch in range(resume_epoch + 1, epochs + 1):
        train_losses = train_epoch(
            model, loaders["train"], optimizer, criterion, device, epoch, positive_sample_types,
            presence_positive_theta_min_mrad, presence_low_angle_policy,
            presence_negative_theta_max_mrad,
        )
        validation = predict(
            model, loaders["validation"], criterion, device,
            positive_sample_types=positive_sample_types,
            presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
            presence_low_angle_policy=presence_low_angle_policy,
            presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
        )
        validation_metrics = result_metrics(validation, threshold)
        balanced_regression_value, regression_bin_metrics = balanced_theta_mae_mrad(
            validation, regression_edges
        )
        validation_metrics["balanced_theta_mae_mrad"] = balanced_regression_value
        validation_metrics["theta_bin_metrics"] = regression_bin_metrics
        validation_losses = {"total": validation["loss"]}
        classification_value = float(validation_metrics[classification_monitor])
        regression_value = float(validation_metrics[regression_monitor])
        row = {
            "epoch": epoch, "train": train_losses, "validation": validation_losses,
            "validation_metrics": validation_metrics,
            "checkpoint_monitor": {
                "classification": {"name": classification_monitor, "mode": "max", "value": classification_value},
                "regression": {"name": regression_monitor, "mode": "min", "value": regression_value},
            },
        }
        history.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
        classification_improved = classification_value > best_classification_score + min_delta
        regression_improved = regression_value < best_regression_score - min_delta
        checkpoint_common = {
            "model_state": model.state_dict(), "config": config, "epoch": epoch,
            "validation_loss": float(validation["loss"]),
            "resolved_regression_balance": regression_balance,
        }
        if classification_improved:
            best_classification_score = classification_value
            classification_stale = 0
            classification_checkpoint = {
                **checkpoint_common,
                "checkpoint_monitor": classification_monitor,
                "checkpoint_mode": "max",
                "checkpoint_score": best_classification_score,
            }
            torch.save(classification_checkpoint, classification_checkpoint_path)
            # Backward-compatible production name: score/presence checkpoint.
            torch.save(classification_checkpoint, checkpoint_path)
        else:
            classification_stale += 1
        if regression_improved:
            best_regression_score = regression_value
            regression_stale = 0
            torch.save({
                **checkpoint_common,
                "checkpoint_monitor": regression_monitor,
                "checkpoint_mode": "min",
                "checkpoint_score": best_regression_score,
            }, regression_checkpoint_path)
        else:
            regression_stale += 1
        if classification_stale >= patience and regression_stale >= patience:
            break
    classification_checkpoint = torch.load(
        classification_checkpoint_path, map_location=device, weights_only=False
    )
    regression_checkpoint = torch.load(
        regression_checkpoint_path, map_location=device, weights_only=False
    )
    model.load_state_dict(classification_checkpoint["model_state"])
    classification_validation = predict(
        model, loaders["validation"], criterion, device,
        positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    classification_test = predict(
        model, loaders["test"], criterion, device,
        positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    regression_model = build_model(config).to(device)
    regression_model.load_state_dict(regression_checkpoint["model_state"])
    regression_validation = predict(
        regression_model, loaders["validation"], criterion, device,
        positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    regression_test = predict(
        regression_model, loaders["test"], criterion, device,
        positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    validation = combine_head_predictions(classification_validation, regression_validation)
    test = combine_head_predictions(classification_test, regression_test)
    z_reverse_classification = predict(
        model, loaders["validation"], criterion, device, z_transform="reverse", z_seed=seed,
        positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    z_reverse_regression = predict(
        regression_model, loaders["validation"], criterion, device,
        z_transform="reverse", z_seed=seed,
        positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    z_reverse = combine_head_predictions(z_reverse_classification, z_reverse_regression)
    z_permute_classification = predict(
        model, loaders["validation"], criterion, device, z_transform="permute", z_seed=seed,
        positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    z_permute_regression = predict(
        regression_model, loaders["validation"], criterion, device,
        z_transform="permute", z_seed=seed,
        positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    z_permute = combine_head_predictions(z_permute_classification, z_permute_regression)
    layer_drop_results = {}
    for transform in (
        "drop_one_random", "drop_two_isolated", "drop_two_consecutive", "drop_peak"
    ):
        classification_dropped = predict(
            model,
            loaders["validation"],
            criterion,
            device,
            z_transform=transform,
            z_seed=seed,
            positive_sample_types=positive_sample_types,
            presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
            presence_low_angle_policy=presence_low_angle_policy,
            presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
        )
        regression_dropped = predict(
            regression_model,
            loaders["validation"],
            criterion,
            device,
            z_transform=transform,
            z_seed=seed,
            positive_sample_types=positive_sample_types,
            presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
            presence_low_angle_policy=presence_low_angle_policy,
            presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
        )
        layer_drop_results[transform] = combine_head_predictions(
            classification_dropped, regression_dropped
        )
    summary = {
        "mode": mode,
        "task": {
            "label_source": "hard_negative_and_signal_theta",
            "include_sample_types": list(include_sample_types),
            "positive_sample_types": list(positive_sample_types),
            "excluded_sample_types": [value for value in SAMPLE_TYPE_NAMES if value not in include_sample_types],
            "signal_theta_range_mrad": None if signal_theta_range_mrad is None else list(signal_theta_range_mrad),
            "presence_positive_theta_min_mrad": presence_positive_theta_min_mrad,
            "presence_low_angle_policy": presence_low_angle_policy,
            "presence_negative_theta_max_mrad": presence_negative_theta_max_mrad,
            "signal_below_presence_minimum": (
                (
                    "bce_negative_below_presence_negative_theta_max_and_"
                    "ignored_until_presence_positive_theta_min_used_for_regression"
                )
                if presence_negative_theta_max_mrad is not None
                else "bce_negative_and_used_for_regression"
                if presence_low_angle_policy == "negative"
                else "ignored_for_bce_used_for_regression"
            ),
        },
        "device": str(device),
        "parameter_count": model.parameter_count,
        "elapsed_seconds": time.time() - started,
        "subset_sizes": {key: len(value.dataset) for key, value in loaders.items()},
        "best_epoch": int(classification_checkpoint["epoch"]),
        "best_regression_epoch": int(regression_checkpoint["epoch"]),
        "checkpoint_selection": {
            "classification": {
                "path": str(classification_checkpoint_path),
                "monitor": classification_checkpoint["checkpoint_monitor"],
                "mode": "max",
                "best_score": float(classification_checkpoint["checkpoint_score"]),
                "epoch": int(classification_checkpoint["epoch"]),
            },
            "regression": {
                "path": str(regression_checkpoint_path),
                "monitor": regression_checkpoint["checkpoint_monitor"],
                "mode": "min",
                "best_score": float(regression_checkpoint["checkpoint_score"]),
                "epoch": int(regression_checkpoint["epoch"]),
            },
        },
        "regression_balance": regression_balance,
        "train_sample_type_fractions": train_sample_type_fractions,
        "history": history,
        "validation": result_metrics(validation, threshold),
        "test": result_metrics(test, threshold),
        "classification_checkpoint_validation": result_metrics(classification_validation, threshold),
        "classification_checkpoint_test": result_metrics(classification_test, threshold),
        "regression_checkpoint_validation": result_metrics(regression_validation, threshold),
        "regression_checkpoint_test": result_metrics(regression_test, threshold),
        "validation_groups": grouped_diagnostics(validation, threshold),
        "z_order_diagnostic": {
            "original": result_metrics(validation, threshold),
            "reversed": result_metrics(z_reverse, threshold),
            "permuted": result_metrics(z_permute, threshold),
            "mean_probability_abs_change_reversed": float(np.abs(validation["probabilities"] - z_reverse["probabilities"]).mean()),
            "mean_probability_abs_change_permuted": float(np.abs(validation["probabilities"] - z_permute["probabilities"]).mean()),
            "mean_slope_l2_change_reversed": float(np.linalg.norm(validation["slope_prediction"] - z_reverse["slope_prediction"], axis=1).mean()),
            "mean_slope_l2_change_permuted": float(np.linalg.norm(validation["slope_prediction"] - z_permute["slope_prediction"], axis=1).mean()),
        },
        "layer_drop_diagnostic": {
            transform: {
                "metrics": result_metrics(result, threshold),
                "mean_probability_abs_change": float(np.abs(
                    validation["probabilities"] - result["probabilities"]
                ).mean()),
                "mean_slope_l2_change": float(np.linalg.norm(
                    validation["slope_prediction"] - result["slope_prediction"], axis=1
                ).mean()),
            }
            for transform, result in layer_drop_results.items()
        },
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
    save_plots(output_dir, history, validation)
    return summary


def run_tiny_overfit(config: dict[str, Any]) -> dict[str, Any]:
    tiny = json.loads(json.dumps(config))
    tiny["subsets"]["tiny"]["validation"] = tiny["subsets"]["tiny"]["train"]
    tiny["subsets"]["tiny"]["test"] = tiny["subsets"]["tiny"]["train"]
    # Fixed central crops make this a genuine memorization test.
    tiny["data"]["dataset"].update({
        "max_jitter": 0, "positive_max_jitter": 0,
        "random_rotations": False, "reflect_x": False, "reflect_y": False,
        "layer_drop_one_probability": 0.0,
        "layer_drop_two_consecutive_probability": 0.0,
    })
    # Use the train split for all three loaders by invoking a compact dedicated loop.
    seed = int(tiny["seed"]); set_seed(seed); device = choose_device(str(tiny.get("device", "auto")))
    path = Path(tiny["data"]["hdf5_path"]).resolve(); output_dir = Path(tiny["output_dir"]).resolve() / "tiny"; output_dir.mkdir(parents=True, exist_ok=True)
    data_cfg = CNNDataConfig(**tiny["data"]["dataset"], seed=seed)
    include_sample_types = tuple(int(value) for value in tiny["task"]["include_sample_types"])
    positive_sample_types = tuple(int(value) for value in tiny["task"]["positive_sample_types"])
    theta_values = tiny["task"].get("signal_theta_range_mrad")
    signal_theta_range_mrad = None if theta_values is None else (float(theta_values[0]), float(theta_values[1]))
    presence_theta_value = tiny["task"].get("presence_positive_theta_min_mrad")
    presence_positive_theta_min_mrad = None if presence_theta_value is None else float(presence_theta_value)
    presence_low_angle_policy = str(
        tiny["task"].get("presence_low_angle_policy", "ignore")
    )
    negative_theta_value = tiny["task"].get("presence_negative_theta_max_mrad")
    presence_negative_theta_max_mrad = (
        None if negative_theta_value is None else float(negative_theta_value)
    )
    size = int(tiny["subsets"]["tiny"]["train"])
    tiny_dataset = HDF5CNN21DDataset(path, split="train", training=False, config=data_cfg)
    with h5py.File(path, "r") as hdf5:
        tiny_sample_types = np.asarray(hdf5["sample_type"][tiny_dataset.indices])
        tiny_slope_xy = np.asarray(hdf5["slope_xy"][tiny_dataset.indices], dtype=np.float64)
    tiny_signal = np.isin(tiny_sample_types, positive_sample_types)
    tiny_theta_mrad = 1000.0 * np.arctan(np.linalg.norm(tiny_slope_xy, axis=1) / 27.0)
    tiny_labels = tiny_signal.copy()
    if presence_positive_theta_min_mrad is not None:
        tiny_labels &= tiny_theta_mrad > presence_positive_theta_min_mrad
    eligible = np.isin(tiny_sample_types, include_sample_types)
    if signal_theta_range_mrad is not None:
        eligible &= ~tiny_signal | (
            (tiny_theta_mrad > signal_theta_range_mrad[0]) & (tiny_theta_mrad < signal_theta_range_mrad[1])
        )
    # Legacy policy excludes low-angle signal. A negative lower band keeps only
    # its explicitly negative portion in the memorization sample.
    if (
        presence_positive_theta_min_mrad is not None
        and presence_low_angle_policy == "ignore"
    ):
        eligible &= ~tiny_signal | tiny_labels
    if (
        presence_low_angle_policy == "negative"
        and presence_negative_theta_max_mrad is not None
    ):
        low_angle_negative = tiny_signal & (
            tiny_theta_mrad < presence_negative_theta_max_mrad
        )
        ignored_transition = tiny_signal & ~tiny_labels & ~low_angle_negative
        eligible &= ~ignored_transition
    tiny_rng = np.random.default_rng(seed)
    if presence_low_angle_policy == "negative" and presence_positive_theta_min_mrad is not None:
        low_angle_signal = np.flatnonzero(eligible & tiny_signal & ~tiny_labels)
        other_negative = np.flatnonzero(eligible & ~tiny_signal)
        low_count = size // 4
        negative = np.concatenate((
            tiny_rng.choice(low_angle_signal, size=low_count, replace=False),
            tiny_rng.choice(other_negative, size=size // 2 - low_count, replace=False),
        ))
    else:
        negative = tiny_rng.choice(
            np.flatnonzero(eligible & ~tiny_labels), size=size // 2, replace=False
        )
    positive = tiny_rng.choice(np.flatnonzero(eligible & tiny_labels), size=size - size // 2, replace=False)
    tiny_indices = np.concatenate((negative, positive)); tiny_rng.shuffle(tiny_indices)
    loader = DataLoader(Subset(tiny_dataset, tiny_indices.tolist()), batch_size=size, shuffle=False, num_workers=0)
    model = build_model(tiny).to(device); criterion = MultiTaskLoss(**tiny["loss"]).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=float(tiny["subsets"]["tiny"]["learning_rate"]), weight_decay=0.0)
    threshold = float(tiny["evaluation"]["threshold"]); history = []
    initial = predict(
        model, loader, criterion, device, positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    epochs = int(tiny["subsets"]["tiny"]["epochs"])
    for epoch in range(1, epochs + 1):
        losses = train_epoch(
            model, loader, optimizer, criterion, device, epoch, positive_sample_types,
            presence_positive_theta_min_mrad, presence_low_angle_policy,
            presence_negative_theta_max_mrad,
        )
        result = predict(
            model, loader, criterion, device, positive_sample_types=positive_sample_types,
            presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
            presence_low_angle_policy=presence_low_angle_policy,
            presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
        )
        metrics = result_metrics(result, threshold)
        history.append({"epoch": epoch, "train": losses, "metrics": metrics})
        if epoch == 1 or epoch % 5 == 0 or epoch == epochs:
            print(json.dumps(history[-1], sort_keys=True), flush=True)
    final = predict(
        model, loader, criterion, device, positive_sample_types=positive_sample_types,
        presence_positive_theta_min_mrad=presence_positive_theta_min_mrad,
        presence_low_angle_policy=presence_low_angle_policy,
        presence_negative_theta_max_mrad=presence_negative_theta_max_mrad,
    )
    summary = {
        "device": str(device), "parameter_count": model.parameter_count, "examples": size,
        "task": {"label_source": "hard_negative_and_signal_theta", "include_sample_types": list(include_sample_types), "positive_sample_types": list(positive_sample_types), "signal_theta_range_mrad": None if signal_theta_range_mrad is None else list(signal_theta_range_mrad), "presence_positive_theta_min_mrad": presence_positive_theta_min_mrad, "presence_low_angle_policy": presence_low_angle_policy, "presence_negative_theta_max_mrad": presence_negative_theta_max_mrad},
        "initial_loss": initial["loss"], "final_loss": final["loss"],
        "initial_metrics": result_metrics(initial, threshold),
        "final_metrics": result_metrics(final, threshold),
        "history": history,
    }
    torch.save({"model_state": model.state_dict(), "config": tiny}, output_dir / "overfit_model.pt")
    with (output_dir / "summary.json").open("w", encoding="utf-8") as stream: json.dump(summary, stream, indent=2)
    # Adapt history for the shared plotter.
    plot_history = [{"epoch": row["epoch"], "train": row["train"], "validation": {"total": row["train"]["total"]}, "validation_metrics": row["metrics"]} for row in history]
    save_plots(output_dir, plot_history, final)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("training/configs/cnn21d.yaml"))
    parser.add_argument("--mode", choices=("tiny", "pilot", "intermediate", "full"), required=True)
    args = parser.parse_args(); config = load_config(args.config)
    summary = run_tiny_overfit(config) if args.mode == "tiny" else run_training(config, args.mode)
    print(json.dumps({key: summary[key] for key in summary if key not in {"history", "validation_groups", "z_order_diagnostic"}}, indent=2))


if __name__ == "__main__":
    main()
