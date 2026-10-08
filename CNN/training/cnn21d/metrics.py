from __future__ import annotations

import math
from typing import Any

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    precision_recall_fscore_support,
    roc_auc_score,
    roc_curve,
)


NEGATIVE_TYPE_NAMES = {-1: "not_applicable", 0: "poisson", 1: "hard", 2: "unknown"}
SAMPLE_TYPE_NAMES = {0: "poisson", 1: "hard", 2: "signal"}


def slopes_to_angles(slope_xy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    slope_xy = np.asarray(slope_xy, dtype=np.float64)
    theta_mrad = 1000.0 * np.arctan(np.linalg.norm(slope_xy, axis=-1) / 27.0)
    phi_rad = np.arctan2(slope_xy[..., 1], slope_xy[..., 0])
    return theta_mrad, phi_rad


def _finite_or_none(value: float) -> float | None:
    return float(value) if np.isfinite(value) else None


def classification_metrics(labels: np.ndarray, probabilities: np.ndarray, threshold: float) -> dict[str, Any]:
    labels = np.asarray(labels, dtype=np.int64)
    probabilities = np.asarray(probabilities, dtype=np.float64)
    predictions = (probabilities >= threshold).astype(np.int64)
    precision, recall, f1, _ = precision_recall_fscore_support(
        labels, predictions, average="binary", zero_division=0
    )
    matrix = confusion_matrix(labels, predictions, labels=[0, 1])
    result: dict[str, Any] = {
        "threshold": float(threshold),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "confusion_matrix": matrix.tolist(),
    }
    if np.unique(labels).size == 2:
        result["auroc"] = float(roc_auc_score(labels, probabilities))
        result["auprc"] = float(average_precision_score(labels, probabilities))
    else:
        result["auroc"] = None
        result["auprc"] = None
    return result


def regression_metrics(labels: np.ndarray, truth: np.ndarray, prediction: np.ndarray) -> dict[str, Any]:
    positive = np.asarray(labels) > 0.5
    if not np.any(positive):
        return {"positive_count": 0}
    truth = np.asarray(truth)[positive]
    prediction = np.asarray(prediction)[positive]
    residual = prediction - truth
    theta_true, phi_true = slopes_to_angles(truth)
    theta_pred, phi_pred = slopes_to_angles(prediction)
    theta_residual = theta_pred - theta_true
    phi_residual = np.arctan2(np.sin(phi_pred - phi_true), np.cos(phi_pred - phi_true))
    return {
        "positive_count": int(positive.sum()),
        "slope_mae": float(np.abs(residual).mean()),
        "slope_mae_sx": float(np.abs(residual[:, 0]).mean()),
        "slope_mae_sy": float(np.abs(residual[:, 1]).mean()),
        "slope_bias_sx": float(residual[:, 0].mean()),
        "slope_bias_sy": float(residual[:, 1].mean()),
        "theta_bias_mrad": float(theta_residual.mean()),
        "theta_mae_mrad": float(np.abs(theta_residual).mean()),
        "theta_resolution_mrad": float(theta_residual.std(ddof=0)),
        "phi_circular_mae_rad": float(np.abs(phi_residual).mean()),
    }


def all_metrics(
    labels: np.ndarray,
    probabilities: np.ndarray,
    slope_truth: np.ndarray,
    slope_prediction: np.ndarray,
    threshold: float,
) -> dict[str, Any]:
    return {
        **classification_metrics(labels, probabilities, threshold),
        **regression_metrics(labels, slope_truth, slope_prediction),
    }


def curve_points(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, np.ndarray]:
    if np.unique(labels).size < 2:
        return {}
    fpr, tpr, roc_thresholds = roc_curve(labels, probabilities)
    precision, recall, pr_thresholds = precision_recall_curve(labels, probabilities)
    return {
        "fpr": fpr,
        "tpr": tpr,
        "roc_thresholds": roc_thresholds,
        "pr_precision": precision,
        "pr_recall": recall,
        "pr_thresholds": pr_thresholds,
    }
