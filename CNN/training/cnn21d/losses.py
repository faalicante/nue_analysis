from __future__ import annotations

import torch
from torch import nn
from torch.nn import functional as F


class MultiTaskLoss(nn.Module):
    """Masked BCE and masked Smooth-L1, each normalized by its effective count."""

    def __init__(
        self,
        *,
        presence_weight: float = 1.0,
        slope_weight: float = 1.0,
        smooth_l1_beta: float = 1.0,
        positive_class_weight: float | None = None,
        regression_theta_bin_edges_mrad: list[float] | tuple[float, ...] | None = None,
        regression_theta_bin_weights: list[float] | tuple[float, ...] | None = None,
    ) -> None:
        super().__init__()
        if presence_weight < 0 or slope_weight < 0:
            raise ValueError("loss weights must be non-negative")
        self.presence_weight = float(presence_weight)
        self.slope_weight = float(slope_weight)
        self.smooth_l1_beta = float(smooth_l1_beta)
        value = [] if positive_class_weight is None else [float(positive_class_weight)]
        self.register_buffer("positive_class_weight", torch.tensor(value, dtype=torch.float32))
        edges = [] if regression_theta_bin_edges_mrad is None else [
            float(value) for value in regression_theta_bin_edges_mrad
        ]
        weights = [] if regression_theta_bin_weights is None else [
            float(value) for value in regression_theta_bin_weights
        ]
        if bool(edges) != bool(weights):
            raise ValueError("regression theta bin edges and weights must be provided together")
        if edges:
            if len(edges) < 2 or any(b <= a for a, b in zip(edges, edges[1:])):
                raise ValueError("regression theta bin edges must be strictly increasing")
            if len(weights) != len(edges) - 1:
                raise ValueError("regression theta bin weights must have len(edges)-1 values")
            if any(value <= 0 for value in weights):
                raise ValueError("regression theta bin weights must be positive")
        self.register_buffer("regression_theta_bin_edges_mrad", torch.tensor(edges))
        self.register_buffer("regression_theta_bin_weights", torch.tensor(weights))

    def _regression_weights(self, slope_xy: torch.Tensor) -> torch.Tensor:
        if not self.regression_theta_bin_weights.numel():
            return torch.ones(slope_xy.shape[0], device=slope_xy.device, dtype=slope_xy.dtype)
        theta = 1000.0 * torch.atan(torch.linalg.vector_norm(slope_xy, dim=-1) / 27.0)
        internal_edges = self.regression_theta_bin_edges_mrad[1:-1].to(
            device=theta.device, dtype=theta.dtype
        )
        bins = torch.bucketize(theta, internal_edges)
        return self.regression_theta_bin_weights.to(
            device=theta.device, dtype=slope_xy.dtype
        )[bins]

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        presence: torch.Tensor,
        slope_xy: torch.Tensor,
        *,
        presence_mask: torch.Tensor | None = None,
        regression_mask: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        if presence_mask is None:
            presence_mask = torch.ones_like(presence, dtype=torch.bool)
        else:
            presence_mask = presence_mask.bool()
        if regression_mask is None:
            regression_mask = presence > 0.5
        else:
            regression_mask = regression_mask.bool()
        pos_weight = self.positive_class_weight if self.positive_class_weight.numel() else None
        labeled_count = presence_mask.sum()
        if bool(presence_mask.any()):
            presence_loss = F.binary_cross_entropy_with_logits(
                outputs["presence_logit"][presence_mask], presence.float()[presence_mask],
                pos_weight=pos_weight,
            )
        else:
            presence_loss = outputs["presence_logit"].sum() * 0.0
        regression_count = regression_mask.sum()
        if bool(regression_mask.any()):
            # Average components per example, then use the configured theta-bin
            # weights and normalize by their effective sum. With unit weights
            # this is exactly normalization by the number of positive examples.
            component_loss = F.smooth_l1_loss(
                outputs["slope_xy"][regression_mask],
                slope_xy[regression_mask],
                reduction="none",
                beta=self.smooth_l1_beta,
            ).mean(dim=-1)
            regression_weights = self._regression_weights(slope_xy[regression_mask])
            slope_loss = (component_loss * regression_weights).sum() / regression_weights.sum()
            effective_regression_weight = regression_weights.sum()
        else:
            # Keeps the graph connected and gives exactly zero slope-head gradients.
            slope_loss = outputs["slope_xy"].sum() * 0.0
            effective_regression_weight = slope_loss.detach()
        total = self.presence_weight * presence_loss + self.slope_weight * slope_loss
        return {
            "total": total,
            "presence": presence_loss,
            "slope": slope_loss,
            "labeled_count": labeled_count,
            "positive_count": regression_count,
            "effective_regression_weight": effective_regression_weight,
        }
