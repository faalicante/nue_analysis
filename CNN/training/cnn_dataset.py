"""PyTorch input pipeline for the 2+1D CNN HDF5 files.

Axis convention
---------------
``crop_bw.C`` reads ROOT histograms as ``GetBinContent(sourceX, sourceY)`` and
writes the flat C array at ``(z * ny * nx) + (y * nx) + x``.  The converter
preserves that memory layout as ``regions_raw[event, z, y, x]``.  Therefore a
NumPy plane is indexed ``plane[y, x]``: columns are ROOT x and rows are ROOT y.

For an image shown with ``origin="lower"``, increasing NumPy row index is
increasing physical y.  A physical counter-clockwise quarter turn is consequently
``np.rot90(..., k=-1)`` (not ``k=+1``).  It maps both an image displacement and a
slope as ``(x, y) -> (-y, x)``.  ``reflect_x`` mirrors the x coordinate/columns
and ``reflect_y`` mirrors the y coordinate/rows.
"""

from __future__ import annotations

import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import h5py
import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset, get_worker_info


SplitName = Literal["train", "validation", "test"]
NormalizationMode = Literal[
    "background_mu_sigma",
    "background_mu_residual",
    "mpv_z_sigma",
    "mpv_z_residual",
]
EXPECTED_SOURCE_SHAPE = (57, 32, 32)
DEFAULT_CROP_SIZE = 20
SPLIT_NAMES: tuple[SplitName, ...] = ("train", "validation", "test")


@dataclass(frozen=True)
class CNNDataConfig:
    """Crop, augmentation, normalization, and diagnostic settings.

    ``positive_max_jitter`` is the requested safety control.  If set, it reduces
    translation only for the shower classes (``sample_type`` 1 or 2); it never
    changes the stored target.  ``min_positive_visible_fraction`` defines the warning flag as
    retained positive excess relative to the central crop, without a sigma cut.
    """

    crop_size: int = DEFAULT_CROP_SIZE
    max_jitter: int = 5
    positive_max_jitter: int | None = None
    random_rotations: bool = True
    reflect_x: bool = True
    reflect_y: bool = True
    normalization_epsilon: float = 1.0e-6
    normalization_mode: NormalizationMode = "background_mu_sigma"
    min_positive_visible_fraction: float = 0.70
    layer_drop_one_probability: float = 0.0
    layer_drop_two_consecutive_probability: float = 0.0
    poisson_background_augmentation_probability: float = 0.0
    poisson_background_mu_ranges: tuple[tuple[float, float], ...] = ()
    poisson_background_mu_weights: tuple[float, ...] = ()
    poisson_background_evaluation_target_mu: float | None = None
    seed: int = 12345

    def __post_init__(self) -> None:
        if self.crop_size <= 0:
            raise ValueError("crop_size must be positive")
        maximum = (EXPECTED_SOURCE_SHAPE[-1] - self.crop_size) // 2
        if self.crop_size > EXPECTED_SOURCE_SHAPE[-1]:
            raise ValueError("crop_size cannot exceed the 32x32 source region")
        if not 0 <= self.max_jitter <= maximum:
            raise ValueError(f"max_jitter must be in [0, {maximum}]")
        if self.positive_max_jitter is not None and not (
            0 <= self.positive_max_jitter <= self.max_jitter
        ):
            raise ValueError("positive_max_jitter must be in [0, max_jitter]")
        if self.normalization_epsilon <= 0:
            raise ValueError("normalization_epsilon must be positive")
        if self.normalization_mode not in {
            "background_mu_sigma",
            "background_mu_residual",
            "mpv_z_sigma",
            "mpv_z_residual",
        }:
            raise ValueError(f"unsupported normalization_mode {self.normalization_mode!r}")
        if not 0.0 <= self.min_positive_visible_fraction <= 1.0:
            raise ValueError("min_positive_visible_fraction must be in [0, 1]")
        if not 0.0 <= self.layer_drop_one_probability <= 1.0:
            raise ValueError("layer_drop_one_probability must be in [0, 1]")
        if not 0.0 <= self.layer_drop_two_consecutive_probability <= 1.0:
            raise ValueError("layer_drop_two_consecutive_probability must be in [0, 1]")
        if (
            self.layer_drop_one_probability
            + self.layer_drop_two_consecutive_probability
            > 1.0
        ):
            raise ValueError("layer-drop probabilities must sum to at most one")
        if not 0.0 <= self.poisson_background_augmentation_probability <= 1.0:
            raise ValueError(
                "poisson_background_augmentation_probability must be in [0, 1]"
            )
        if self.poisson_background_augmentation_probability > 0.0:
            if self.normalization_mode.startswith("mpv_z"):
                raise ValueError(
                    "Poisson background augmentation is unsupported with local MPV_z normalization"
                )
            if not self.poisson_background_mu_ranges:
                raise ValueError(
                    "poisson_background_mu_ranges is required when Poisson augmentation is enabled"
                )
            for bounds in self.poisson_background_mu_ranges:
                if len(bounds) != 2:
                    raise ValueError("each Poisson background mu range needs [minimum, maximum]")
                minimum, maximum = (float(bounds[0]), float(bounds[1]))
                if minimum <= 0.0 or maximum < minimum:
                    raise ValueError("invalid Poisson background mu range")
            if self.poisson_background_mu_weights:
                if len(self.poisson_background_mu_weights) != len(
                    self.poisson_background_mu_ranges
                ):
                    raise ValueError(
                        "poisson_background_mu_weights must match the number of ranges"
                    )
                if any(float(weight) < 0.0 for weight in self.poisson_background_mu_weights):
                    raise ValueError("Poisson background mu weights must be non-negative")
                if sum(float(weight) for weight in self.poisson_background_mu_weights) <= 0.0:
                    raise ValueError("Poisson background mu weights must have a positive sum")
        if (
            self.poisson_background_evaluation_target_mu is not None
            and self.poisson_background_evaluation_target_mu <= 0.0
        ):
            raise ValueError("poisson_background_evaluation_target_mu must be positive")
        if self.seed < 0:
            raise ValueError("seed must be non-negative")

    def jitter_for(self, is_shower_class: bool) -> int:
        if is_shower_class and self.positive_max_jitter is not None:
            return self.positive_max_jitter
        return self.max_jitter


@dataclass(frozen=True)
class CNNLoaderConfig:
    batch_size: int = 32
    num_workers: int = 0
    pin_memory: bool = False
    persistent_workers: bool = False
    prefetch_factor: int | None = None
    shuffle_train: bool = True
    drop_last_train: bool = False
    seed: int = 12345

    def __post_init__(self) -> None:
        if self.batch_size <= 0:
            raise ValueError("batch_size must be positive")
        if self.num_workers < 0:
            raise ValueError("num_workers must be non-negative")
        if self.persistent_workers and self.num_workers == 0:
            raise ValueError("persistent_workers requires num_workers > 0")
        if self.prefetch_factor is not None and self.prefetch_factor <= 0:
            raise ValueError("prefetch_factor must be positive")
        if self.seed < 0:
            raise ValueError("seed must be non-negative")


def extract_xy_crop(
    volume: np.ndarray,
    *,
    crop_size: int = DEFAULT_CROP_SIZE,
    offset_xy: tuple[int, int] = (0, 0),
) -> np.ndarray:
    """Return a regular slice crop; out-of-range offsets fail instead of wrapping."""

    if volume.ndim != 3:
        raise ValueError(f"volume must have [z,y,x] shape, found {volume.shape}")
    _, height, width = volume.shape
    if crop_size <= 0 or crop_size > min(height, width):
        raise ValueError("invalid crop_size")
    offset_x, offset_y = (int(offset_xy[0]), int(offset_xy[1]))
    first_x = (width - crop_size) // 2 + offset_x
    first_y = (height - crop_size) // 2 + offset_y
    last_x = first_x + crop_size
    last_y = first_y + crop_size
    if first_x < 0 or first_y < 0 or last_x > width or last_y > height:
        raise ValueError(
            f"crop offset {offset_xy} puts [{first_y}:{last_y},{first_x}:{last_x}] "
            f"outside source shape {(height, width)}"
        )
    return volume[:, first_y:last_y, first_x:last_x]


def normalize_per_slice(
    crop: np.ndarray,
    background_mu: float | np.floating[Any],
    *,
    epsilon: float = 1.0e-6,
    mode: NormalizationMode = "background_mu_sigma",
    mpv_z: np.ndarray | None = None,
    sigma_z: np.ndarray | None = None,
) -> np.ndarray:
    """Normalize raw crop counts using a scalar or local per-z reference."""

    values = np.asarray(crop, dtype=np.float32)
    mu_array = np.asarray(background_mu, dtype=np.float32)
    if values.ndim != 3 or mu_array.shape != ():
        raise ValueError(
            "expected crop [z,y,x] and scalar background_mu, found "
            f"{values.shape} and {mu_array.shape}"
        )
    mu = np.float32(mu_array)
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError("background_mu must be finite and positive")
    if mode == "background_mu_sigma":
        denominator = np.float32(np.sqrt(mu + np.float32(epsilon)))
        return (values - mu) / denominator
    if mode == "background_mu_residual":
        return values - mu
    if mode not in {"mpv_z_sigma", "mpv_z_residual"}:
        raise ValueError(f"unsupported normalization mode {mode!r}")
    if mpv_z is None:
        raise ValueError(f"{mode} requires mpv_z")
    local_mpv = np.asarray(mpv_z, dtype=np.float32)
    if local_mpv.shape != (values.shape[0],) or not np.isfinite(local_mpv).all():
        raise ValueError("mpv_z must be finite with one value per z slice")
    residual = values - local_mpv[:, None, None]
    if mode == "mpv_z_residual":
        return residual
    if sigma_z is None:
        raise ValueError("mpv_z_sigma requires sigma_z")
    local_sigma = np.asarray(sigma_z, dtype=np.float32)
    if (
        local_sigma.shape != (values.shape[0],)
        or not np.isfinite(local_sigma).all()
        or np.any(local_sigma <= 0.0)
    ):
        raise ValueError("sigma_z must be finite and positive with one value per z slice")
    return residual / (local_sigma[:, None, None] + np.float32(epsilon))


def local_mpv_sigma_per_z(volume: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Estimate local background mode and width for every z layer.

    The modal count is insensitive to sparse shower pixels.  The width is the
    RMS of the lower side around that mode, so high-count shower tails cannot
    inflate it.  It is measured on the unjittered 32x32 source region.
    """

    values = np.asarray(volume)
    if values.ndim != 3:
        raise ValueError(f"expected volume [z,y,x], found {values.shape}")
    if values.dtype.kind not in "iu" or np.any(values < 0):
        raise ValueError("local MPV statistics require non-negative integer counts")
    mpv = np.empty(values.shape[0], dtype=np.float32)
    sigma = np.empty(values.shape[0], dtype=np.float32)
    for z, layer in enumerate(values):
        flattened = layer.reshape(-1)
        mode = int(np.argmax(np.bincount(flattened)))
        lower_side = flattened[flattened <= mode].astype(np.float32) - np.float32(mode)
        width = float(np.sqrt(np.mean(np.square(lower_side), dtype=np.float64)))
        # A perfectly constant layer contains no measurable lower tail.  Keep a
        # finite scale rather than allowing a single layer to dominate a crop.
        if width <= 0.0:
            width = float(np.sqrt(max(mode, 1)))
        mpv[z] = mode
        sigma[z] = width
    return mpv, sigma


def transform_xy(
    volume: np.ndarray,
    slope_xy: np.ndarray,
    *,
    quarter_turns_ccw: int = 0,
    reflect_x: bool = False,
    reflect_y: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform an image and ``(slope_x, slope_y)`` in physical ROOT xy.

    Rotation is applied first, then x reflection, then y reflection.  Reflections
    name the coordinate whose sign changes: x reflection reverses columns and
    y reflection reverses rows.
    """

    if volume.ndim < 2:
        raise ValueError("volume must have at least two dimensions")
    slope = np.asarray(slope_xy, dtype=np.float32).reshape(2).copy()
    turns = int(quarter_turns_ccw) % 4

    # With [y,x] indexing and origin="lower", k=-1 is physical xy CCW.
    transformed = np.rot90(volume, k=-turns, axes=(-2, -1))
    for _ in range(turns):
        slope[0], slope[1] = -slope[1], slope[0]
    if reflect_x:
        transformed = transformed[..., ::-1]
        slope[0] = -slope[0]
    if reflect_y:
        transformed = transformed[..., ::-1, :]
        slope[1] = -slope[1]
    return np.ascontiguousarray(transformed), np.ascontiguousarray(slope)


def _decode_hdf5_string(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


class HDF5CNN21DDataset(Dataset[dict[str, Any]]):
    """Lazy, worker-safe dataset over one logical HDF5 split."""

    def __init__(
        self,
        path: str | Path,
        *,
        split: SplitName,
        training: bool | None = None,
        config: CNNDataConfig | None = None,
    ) -> None:
        self.path = Path(path).expanduser().resolve()
        self.split = split
        if split not in SPLIT_NAMES:
            raise ValueError(f"split must be one of {SPLIT_NAMES}")
        self.training = split == "train" if training is None else bool(training)
        if self.training and split != "train":
            raise ValueError("random augmentation is allowed only for the training split")
        self.config = config or CNNDataConfig()
        self._epoch = 0
        self._hdf5: h5py.File | None = None
        self._handle_identity: tuple[int, int] | None = None

        # Schema/index inspection is eager but this handle is closed immediately.
        # Actual samples are opened lazily in the process/worker that reads them.
        with h5py.File(self.path, "r") as hdf5:
            self._validate_schema(hdf5)
            self.indices = np.asarray(hdf5[f"split_indices/{split}"][:], dtype=np.int64)
            if self.indices.size:
                split_code = SPLIT_NAMES.index(split)
                if np.any(np.asarray(hdf5["split"][self.indices]) != split_code):
                    raise ValueError(f"split_indices/{split} disagrees with the split dataset")

    @staticmethod
    def _validate_schema(hdf5: h5py.File) -> None:
        required = {
            "regions_raw",
            "background_mu",
            "presence",
            "slope_xy",
            "sample_type",
            "signal_event_id",
            "tile_id",
            "crop_id_in_tile",
            "negative_type",
            "split",
            "source_file",
            "source_entry",
            "split_indices",
        }
        missing = sorted(required.difference(hdf5.keys()))
        if missing:
            raise ValueError(f"missing HDF5 datasets/groups: {missing}")
        if hdf5["regions_raw"].shape[1:] != EXPECTED_SOURCE_SHAPE:
            raise ValueError(
                f"regions_raw must have [N,57,32,32], found {hdf5['regions_raw'].shape}"
            )
        count = hdf5["regions_raw"].shape[0]
        if hdf5["background_mu"].shape not in {(57,), (count, 57)}:
            raise ValueError(
                "background_mu must have [57] or [N,57], found "
                f"{hdf5['background_mu'].shape}"
            )
        background = np.asarray(hdf5["background_mu"][:], dtype=np.float32)
        if background.ndim == 1:
            varies_with_z = np.any(background != background[0])
        else:
            varies_with_z = np.any(background != background[:, :1])
        if varies_with_z:
            raise ValueError(
                "background_mu is physically scalar per event but the HDF5 values vary with z"
            )
        if hdf5["slope_xy"].shape != (count, 2):
            raise ValueError(f"slope_xy must have [N,2], found {hdf5['slope_xy'].shape}")
        for split in SPLIT_NAMES:
            if split not in hdf5["split_indices"]:
                raise ValueError(f"missing split_indices/{split}")

    def __len__(self) -> int:
        return int(self.indices.size)

    def set_epoch(self, epoch: int) -> None:
        """Select a deterministic, different augmentation stream for an epoch."""

        if epoch < 0:
            raise ValueError("epoch must be non-negative")
        self._epoch = int(epoch)

    def _identity(self) -> tuple[int, int]:
        worker = get_worker_info()
        return os.getpid(), -1 if worker is None else worker.id

    def _file(self) -> h5py.File:
        identity = self._identity()
        if self._hdf5 is None or self._handle_identity != identity:
            self.close()
            self._hdf5 = h5py.File(self.path, "r")
            self._handle_identity = identity
        return self._hdf5

    def close(self) -> None:
        if self._hdf5 is not None:
            try:
                self._hdf5.close()
            finally:
                self._hdf5 = None
                self._handle_identity = None

    def __getstate__(self) -> dict[str, Any]:
        state = self.__dict__.copy()
        state["_hdf5"] = None
        state["_handle_identity"] = None
        return state

    def __del__(self) -> None:
        self.close()

    def _background_for(self, hdf5: h5py.File, raw_index: int) -> np.float32:
        background = hdf5["background_mu"]
        values = background[:] if background.ndim == 1 else background[raw_index]
        return np.float32(values[0])

    def _augmentation_parameters(
        self, raw_index: int, is_shower_class: bool
    ) -> tuple[int, int, int, bool, bool]:
        if not self.training:
            return 0, 0, 0, False, False
        rng = np.random.default_rng(
            np.random.SeedSequence([self.config.seed, self._epoch, int(raw_index)])
        )
        jitter = self.config.jitter_for(is_shower_class)
        offset_x = int(rng.integers(-jitter, jitter + 1))
        offset_y = int(rng.integers(-jitter, jitter + 1))
        turns = int(rng.integers(0, 4)) if self.config.random_rotations else 0
        flip_x = bool(rng.integers(0, 2)) if self.config.reflect_x else False
        flip_y = bool(rng.integers(0, 2)) if self.config.reflect_y else False
        return offset_x, offset_y, turns, flip_x, flip_y

    def _dropped_layers(self, raw_index: int, z_count: int) -> tuple[int, ...]:
        """Select deterministic train-only missing layers for reconstruction robustness."""

        if not self.training:
            return ()
        rng = np.random.default_rng(
            np.random.SeedSequence([self.config.seed, self._epoch, int(raw_index), 99173])
        )
        draw = float(rng.random())
        if draw < self.config.layer_drop_two_consecutive_probability:
            first = int(rng.integers(0, z_count - 1))
            return first, first + 1
        if draw < (
            self.config.layer_drop_two_consecutive_probability
            + self.config.layer_drop_one_probability
        ):
            return (int(rng.integers(0, z_count)),)
        return ()

    def _augment_poisson_background(
        self,
        crop: np.ndarray,
        background_mu: np.float32,
        raw_index: int,
    ) -> tuple[np.ndarray, np.float32, np.float32]:
        """Add masked Poisson counts to a sampled or fixed evaluation background level."""

        if self.training:
            probability = self.config.poisson_background_augmentation_probability
            if probability <= 0.0:
                return crop, background_mu, np.float32(0.0)
            rng = np.random.default_rng(
                np.random.SeedSequence([self.config.seed, self._epoch, int(raw_index), 74891])
            )
            if float(rng.random()) >= probability:
                return crop, background_mu, np.float32(0.0)
            ranges = self.config.poisson_background_mu_ranges
            weights = self.config.poisson_background_mu_weights
            probabilities = None
            if weights:
                probabilities = np.asarray(weights, dtype=np.float64)
                probabilities /= probabilities.sum()
            range_index = int(rng.choice(len(ranges), p=probabilities))
            minimum, maximum = (float(value) for value in ranges[range_index])
            sampled_target = (
                minimum if minimum == maximum else float(rng.uniform(minimum, maximum))
            )
        else:
            evaluation_target = self.config.poisson_background_evaluation_target_mu
            if evaluation_target is None:
                return crop, background_mu, np.float32(0.0)
            sampled_target = float(evaluation_target)
            rng = np.random.default_rng(
                np.random.SeedSequence([self.config.seed, int(raw_index), 74892])
            )
        target_mu = np.float32(max(float(background_mu), sampled_target))
        added_mu = np.float32(target_mu - background_mu)
        if added_mu <= 0.0:
            return crop, background_mu, np.float32(0.0)
        noise = rng.poisson(float(added_mu), size=crop.shape).astype(np.float32)
        noise[np.asarray(crop) <= 0] = 0.0
        return np.asarray(crop, dtype=np.float32) + noise, target_mu, added_mu

    def __getitem__(self, index: int) -> dict[str, Any]:
        raw_index = int(self.indices[index])
        hdf5 = self._file()
        raw = np.asarray(hdf5["regions_raw"][raw_index])
        original_background_mu = self._background_for(hdf5, raw_index)
        local_mpv_z: np.ndarray | None = None
        local_sigma_z: np.ndarray | None = None
        if self.config.normalization_mode.startswith("mpv_z"):
            cached_mpv = hdf5.get("metadata/local_mpv_z")
            cached_sigma = hdf5.get("metadata/local_sigma_z")
            if cached_mpv is not None and cached_sigma is not None:
                local_mpv_z = np.asarray(cached_mpv[raw_index], dtype=np.float32)
                local_sigma_z = np.asarray(cached_sigma[raw_index], dtype=np.float32)
            else:
                local_mpv_z, local_sigma_z = local_mpv_sigma_per_z(raw)
        presence_value = float(hdf5["presence"][raw_index])
        sample_type_value = int(hdf5["sample_type"][raw_index])
        original_slope = np.asarray(hdf5["slope_xy"][raw_index], dtype=np.float32)
        offset_x, offset_y, turns, flip_x, flip_y = self._augmentation_parameters(
            raw_index, sample_type_value != 0
        )

        # Crop is a plain slice (no roll/padding/wrap), then normalization is
        # performed independently for every z slice, then rigid xy transforms.
        crop = extract_xy_crop(
            raw, crop_size=self.config.crop_size, offset_xy=(offset_x, offset_y)
        )
        crop_for_network, background_mu, added_background_mu = (
            self._augment_poisson_background(crop, original_background_mu, raw_index)
        )
        normalized = normalize_per_slice(
            crop_for_network,
            background_mu,
            epsilon=self.config.normalization_epsilon,
            mode=self.config.normalization_mode,
            mpv_z=local_mpv_z,
            sigma_z=local_sigma_z,
        )
        normalized, slope = transform_xy(
            normalized,
            original_slope,
            quarter_turns_ccw=turns,
            reflect_x=flip_x,
            reflect_y=flip_y,
        )
        dropped_layers = self._dropped_layers(raw_index, normalized.shape[0])
        if dropped_layers:
            # A missing reconstructed layer is represented as zero raw content,
            # expressed in the same normalization used by every other voxel.
            if self.config.normalization_mode == "background_mu_sigma":
                missing_value: float | np.ndarray = np.float32(
                    -background_mu / np.sqrt(background_mu + self.config.normalization_epsilon)
                )
            elif self.config.normalization_mode == "background_mu_residual":
                missing_value = np.float32(-background_mu)
            elif self.config.normalization_mode == "mpv_z_sigma":
                assert local_mpv_z is not None and local_sigma_z is not None
                missing_value = -local_mpv_z / (
                    local_sigma_z + np.float32(self.config.normalization_epsilon)
                )
            else:
                assert local_mpv_z is not None
                missing_value = -local_mpv_z
            if np.ndim(missing_value):
                normalized[list(dropped_layers)] = np.asarray(missing_value)[
                    list(dropped_layers), None, None
                ]
            else:
                normalized[list(dropped_layers)] = missing_value

        visible_fraction = 1.0
        may_be_insufficient = False
        if sample_type_value != 0:
            central = extract_xy_crop(raw, crop_size=self.config.crop_size)
            central_excess = np.maximum(
                central.astype(np.float32) - original_background_mu, 0.0
            ).sum(dtype=np.float64)
            crop_excess = np.maximum(
                crop.astype(np.float32) - original_background_mu, 0.0
            ).sum(dtype=np.float64)
            visible_fraction = float(crop_excess / central_excess) if central_excess > 0 else 0.0
            may_be_insufficient = (
                visible_fraction < self.config.min_positive_visible_fraction
            )

        worker = get_worker_info()
        metadata: dict[str, Any] = {
            "dataset_index": torch.tensor(index, dtype=torch.int64),
            "hdf5_index": torch.tensor(raw_index, dtype=torch.int64),
            "split_code": torch.tensor(int(hdf5["split"][raw_index]), dtype=torch.int64),
            "split_name": self.split,
            "sample_type": torch.tensor(sample_type_value, dtype=torch.int64),
            "presence": torch.tensor(presence_value, dtype=torch.float32),
            "background_mu": torch.tensor(background_mu, dtype=torch.float32),
            "background_mu_original": torch.tensor(
                original_background_mu, dtype=torch.float32
            ),
            "poisson_background_added_mu": torch.tensor(
                added_background_mu, dtype=torch.float32
            ),
            "poisson_background_augmented": torch.tensor(
                added_background_mu > 0.0, dtype=torch.bool
            ),
            "signal_event_id": torch.tensor(
                int(hdf5["signal_event_id"][raw_index]), dtype=torch.int64
            ),
            "tile_id": torch.tensor(int(hdf5["tile_id"][raw_index]), dtype=torch.int64),
            "crop_id_in_tile": torch.tensor(
                int(hdf5["crop_id_in_tile"][raw_index]), dtype=torch.int64
            ),
            "negative_type": torch.tensor(
                int(hdf5["negative_type"][raw_index]), dtype=torch.int64
            ),
            "source_file": _decode_hdf5_string(hdf5["source_file"][raw_index]),
            "source_entry": torch.tensor(
                int(hdf5["source_entry"][raw_index]), dtype=torch.int64
            ),
            "crop_offset_xy": torch.tensor([offset_x, offset_y], dtype=torch.int64),
            "rotation_quarter_turns_ccw": torch.tensor(turns, dtype=torch.int64),
            "reflected_x": torch.tensor(flip_x, dtype=torch.bool),
            "reflected_y": torch.tensor(flip_y, dtype=torch.bool),
            "positive_visible_fraction": torch.tensor(visible_fraction, dtype=torch.float32),
            "positive_may_be_insufficiently_visible": torch.tensor(
                may_be_insufficient, dtype=torch.bool
            ),
            "layer_drop_count": torch.tensor(len(dropped_layers), dtype=torch.int64),
            "layer_drop_first": torch.tensor(
                -1 if not dropped_layers else dropped_layers[0], dtype=torch.int64
            ),
            "worker_id": torch.tensor(-1 if worker is None else worker.id, dtype=torch.int64),
            "worker_pid": torch.tensor(os.getpid(), dtype=torch.int64),
        }
        if "metadata/zScale" in hdf5:
            metadata["z_scale"] = torch.tensor(
                float(hdf5["metadata/zScale"][raw_index]), dtype=torch.float32
            )

        return {
            "volume": torch.from_numpy(normalized[None]).to(dtype=torch.float32),
            "sample_type": torch.tensor(sample_type_value, dtype=torch.int64),
            "slope_xy": torch.from_numpy(slope).to(dtype=torch.float32),
            "metadata": metadata,
        }

    def audit_positive_visibility(
        self,
        *,
        chunk_size: int = 64,
        max_reported_indices: int = 50,
    ) -> dict[str, Any]:
        """Count positives at risk under any configured allowed translation.

        The metric is the worst retained ``max(raw-background_mu, 0)`` mass
        divided by that of the central crop.  It uses no 5-sigma threshold and
        does not edit labels.  Rotations/reflections preserve this mass, so only
        crop translations need to be enumerated.
        """

        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        jitter = self.config.jitter_for(True)
        all_ratios: list[np.ndarray] = []
        at_risk_indices: list[int] = []
        total_positive = 0

        with h5py.File(self.path, "r") as hdf5:
            for start in range(0, len(self), chunk_size):
                raw_indices = self.indices[start : start + chunk_size]
                sample_type = np.asarray(hdf5["sample_type"][raw_indices])
                positive_indices = raw_indices[sample_type != 0]
                if positive_indices.size == 0:
                    continue
                total_positive += int(positive_indices.size)
                raw = np.asarray(hdf5["regions_raw"][positive_indices], dtype=np.float32)
                background = hdf5["background_mu"]
                if background.ndim == 1:
                    mu = np.full(
                        positive_indices.size,
                        np.float32(background[0]),
                        dtype=np.float32,
                    )
                else:
                    mu = np.asarray(background[positive_indices, 0], dtype=np.float32)
                excess_xy = np.maximum(raw - mu[:, None, None, None], 0.0).sum(
                    axis=1, dtype=np.float64
                )

                first = (EXPECTED_SOURCE_SHAPE[-1] - self.config.crop_size) // 2
                central_mass = excess_xy[
                    :, first : first + self.config.crop_size, first : first + self.config.crop_size
                ].sum(axis=(1, 2), dtype=np.float64)
                worst_mass = np.full(positive_indices.size, np.inf, dtype=np.float64)
                for offset_y in range(-jitter, jitter + 1):
                    y0 = first + offset_y
                    for offset_x in range(-jitter, jitter + 1):
                        x0 = first + offset_x
                        mass = excess_xy[
                            :,
                            y0 : y0 + self.config.crop_size,
                            x0 : x0 + self.config.crop_size,
                        ].sum(axis=(1, 2), dtype=np.float64)
                        worst_mass = np.minimum(worst_mass, mass)
                ratios = np.divide(
                    worst_mass,
                    central_mass,
                    out=np.zeros_like(worst_mass),
                    where=central_mass > 0,
                )
                all_ratios.append(ratios)
                risky = positive_indices[
                    ratios < self.config.min_positive_visible_fraction
                ]
                remaining = max_reported_indices - len(at_risk_indices)
                if remaining > 0:
                    at_risk_indices.extend(int(value) for value in risky[:remaining])

        ratios = np.concatenate(all_ratios) if all_ratios else np.empty(0, dtype=np.float64)
        at_risk = int(np.count_nonzero(ratios < self.config.min_positive_visible_fraction))
        quantiles = (
            {str(q): float(np.quantile(ratios, q)) for q in (0.0, 0.05, 0.5, 0.95, 1.0)}
            if ratios.size
            else {}
        )
        return {
            "split": self.split,
            "positive_examples": total_positive,
            "at_risk_examples": at_risk,
            "at_risk_fraction": at_risk / total_positive if total_positive else 0.0,
            "max_jitter_for_positive": jitter,
            "minimum_visible_fraction": self.config.min_positive_visible_fraction,
            "worst_case_visible_fraction_quantiles": quantiles,
            "reported_hdf5_indices": at_risk_indices,
            "sample_type_labels_changed": 0,
            "metric": "worst translated positive-excess mass / central-crop positive-excess mass",
        }


def seed_dataloader_worker(worker_id: int) -> None:
    """Seed libraries used by possible downstream transforms in each worker."""

    del worker_id
    worker_seed = torch.initial_seed() % (2**32)
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def create_cnn_dataloaders(
    path: str | Path,
    *,
    data_config: CNNDataConfig | None = None,
    loader_config: CNNLoaderConfig | None = None,
) -> dict[SplitName, DataLoader[dict[str, Any]]]:
    """Build deterministic train/validation/test loaders from HDF5 split indices."""

    data_config = data_config or CNNDataConfig()
    loader_config = loader_config or CNNLoaderConfig()
    loaders: dict[SplitName, DataLoader[dict[str, Any]]] = {}
    for split_index, split in enumerate(SPLIT_NAMES):
        dataset = HDF5CNN21DDataset(
            path,
            split=split,
            training=split == "train",
            config=data_config,
        )
        generator = torch.Generator()
        generator.manual_seed(loader_config.seed + split_index)
        kwargs: dict[str, Any] = {}
        if loader_config.num_workers > 0 and loader_config.prefetch_factor is not None:
            kwargs["prefetch_factor"] = loader_config.prefetch_factor
        loaders[split] = DataLoader(
            dataset,
            batch_size=loader_config.batch_size,
            shuffle=loader_config.shuffle_train if split == "train" else False,
            num_workers=loader_config.num_workers,
            pin_memory=loader_config.pin_memory,
            persistent_workers=loader_config.persistent_workers,
            drop_last=loader_config.drop_last_train if split == "train" else False,
            worker_init_fn=seed_dataloader_worker,
            generator=generator,
            **kwargs,
        )
    return loaders
