from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import nn


def _valid_groups(channels: int, requested: int) -> int:
    """Return the largest requested-or-smaller divisor of ``channels``."""

    for groups in range(min(channels, requested), 0, -1):
        if channels % groups == 0:
            return groups
    return 1


@dataclass(frozen=True)
class CNN21DConfig:
    in_channels: int = 1
    channels: tuple[int, ...] = (16, 32, 64, 128)
    group_norm_groups: tuple[int, ...] = (4, 8, 8, 16)
    pool_kernels: tuple[tuple[int, int, int], ...] = (
        (1, 2, 2),
        (1, 2, 2),
        (2, 1, 1),
        (2, 1, 1),
    )
    head_hidden: int = 64
    dropout: float = 0.0

    def __post_init__(self) -> None:
        count = len(self.channels)
        if count == 0 or len(self.group_norm_groups) != count or len(self.pool_kernels) != count:
            raise ValueError("channels, group_norm_groups and pool_kernels must have equal non-zero length")
        if any(value <= 0 for value in self.channels + self.group_norm_groups):
            raise ValueError("channel and group counts must be positive")
        if not 0.0 <= self.dropout < 1.0:
            raise ValueError("dropout must be in [0, 1)")


class FactorizedConvBlock(nn.Module):
    """Spatial (1x3x3) convolution followed by a z (3x1x1) convolution."""

    def __init__(self, in_channels: int, out_channels: int, groups: int, pool: tuple[int, int, int]) -> None:
        super().__init__()
        norm_groups = _valid_groups(out_channels, groups)
        self.spatial = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=(1, 3, 3), padding=(0, 1, 1), bias=False),
            nn.GroupNorm(norm_groups, out_channels),
            nn.SiLU(inplace=True),
        )
        self.longitudinal = nn.Sequential(
            nn.Conv3d(out_channels, out_channels, kernel_size=(3, 1, 1), padding=(1, 0, 0), bias=False),
            nn.GroupNorm(norm_groups, out_channels),
            nn.SiLU(inplace=True),
        )
        self.pool = nn.MaxPool3d(kernel_size=pool, stride=pool)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.pool(self.longitudinal(self.spatial(x)))


class TwoHeadCNN21D(nn.Module):
    def __init__(self, config: CNN21DConfig | None = None) -> None:
        super().__init__()
        self.config = config or CNN21DConfig()
        blocks: list[nn.Module] = []
        previous = self.config.in_channels
        for channels, groups, pool in zip(
            self.config.channels, self.config.group_norm_groups, self.config.pool_kernels, strict=True
        ):
            blocks.append(FactorizedConvBlock(previous, channels, groups, pool))
            previous = channels
        self.blocks = nn.ModuleList(blocks)
        feature_size = 2 * previous
        self.presence_head = self._head(feature_size, 1)
        self.slope_head = self._head(feature_size, 2)

    def _head(self, feature_size: int, output_size: int) -> nn.Sequential:
        return nn.Sequential(
            nn.Linear(feature_size, self.config.head_hidden),
            nn.SiLU(inplace=True),
            nn.Dropout(self.config.dropout),
            nn.Linear(self.config.head_hidden, output_size),
        )

    @property
    def parameter_count(self) -> int:
        return sum(parameter.numel() for parameter in self.parameters())

    def extract_features(
        self, x: torch.Tensor, *, return_shapes: bool = False
    ) -> torch.Tensor | tuple[torch.Tensor, dict[str, tuple[int, ...]]]:
        if x.ndim != 5:
            raise ValueError(f"expected [B,C,Z,Y,X], found {tuple(x.shape)}")
        shapes: dict[str, tuple[int, ...]] = {"input": tuple(x.shape)}
        for index, block in enumerate(self.blocks, start=1):
            x = block(x)
            shapes[f"block{index}"] = tuple(x.shape)
        average = x.mean(dim=(2, 3, 4))
        maximum = x.amax(dim=(2, 3, 4))
        features = torch.cat((average, maximum), dim=1)
        shapes["global_features"] = tuple(features.shape)
        if return_shapes:
            return features, shapes
        return features

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        features = self.extract_features(x)
        assert isinstance(features, torch.Tensor)
        return {
            "presence_logit": self.presence_head(features).squeeze(-1),
            "slope_xy": self.slope_head(features),
        }
