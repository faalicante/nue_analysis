from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest
import torch

from training.cnn_dataset import (
    CNNDataConfig,
    CNNLoaderConfig,
    HDF5CNN21DDataset,
    create_cnn_dataloaders,
    extract_xy_crop,
    local_mpv_sigma_per_z,
    normalize_per_slice,
    transform_xy,
)


SOURCE_SHAPE = (57, 32, 32)


def write_hdf5(
    path: Path,
    *,
    regions: np.ndarray | None = None,
    background_mu: np.ndarray | None = None,
    presence: np.ndarray | None = None,
) -> None:
    count = 9
    if regions is None:
        regions = np.full((count, *SOURCE_SHAPE), 5, dtype=np.int16)
        for event in range(count):
            regions[event, :, 8 + event % 3 : 12 + event % 3, 9:13] += 10
    if background_mu is None:
        background_mu = np.full((count, 57), 4.0, dtype=np.float32)
    if presence is None:
        presence = np.asarray([1, 0, 1, 0, 1, 1, 0, 1, 0], dtype=np.uint8)
    split = np.asarray([0, 0, 0, 0, 0, 1, 1, 2, 2], dtype=np.int8)
    sample_type = np.asarray([2, 0, 1, 0, 2, 1, 0, 2, 0], dtype=np.int8)
    slopes = np.column_stack(
        (
            np.linspace(0.1, 0.9, count, dtype=np.float32),
            np.linspace(-0.2, 0.6, count, dtype=np.float32),
        )
    )
    slopes[sample_type != 2] = 0
    string_dtype = h5py.string_dtype(encoding="utf-8")

    with h5py.File(path, "w") as hdf5:
        hdf5.attrs["axis_order"] = "event,z,y,x"
        hdf5.create_dataset("regions_raw", data=regions)
        hdf5.create_dataset("background_mu", data=background_mu)
        hdf5.create_dataset("presence", data=presence)
        hdf5.create_dataset("slope_xy", data=slopes)
        hdf5.create_dataset("sample_type", data=sample_type)
        hdf5.create_dataset("signal_event_id", data=np.arange(count, dtype=np.int64))
        hdf5.create_dataset("tile_id", data=np.full(count, -1, dtype=np.int64))
        hdf5.create_dataset("crop_id_in_tile", data=np.full(count, -1, dtype=np.int64))
        hdf5.create_dataset("negative_type", data=np.full(count, -1, dtype=np.int8))
        hdf5.create_dataset("split", data=split)
        hdf5.create_dataset(
            "source_file",
            data=np.asarray(["synthetic.root"] * count, dtype=object),
            dtype=string_dtype,
        )
        hdf5.create_dataset("source_entry", data=np.arange(count, dtype=np.int64))
        metadata = hdf5.create_group("metadata")
        metadata.create_dataset("zScale", data=np.arange(count, dtype=np.float32))
        splits = hdf5.create_group("split_indices")
        splits.create_dataset("train", data=np.flatnonzero(split == 0))
        splits.create_dataset("validation", data=np.flatnonzero(split == 1))
        splits.create_dataset("test", data=np.flatnonzero(split == 2))


@pytest.fixture
def synthetic_hdf5(tmp_path: Path) -> Path:
    path = tmp_path / "dataset.h5"
    write_hdf5(path)
    return path


def test_shapes_dtypes_and_metadata(synthetic_hdf5: Path) -> None:
    dataset = HDF5CNN21DDataset(synthetic_hdf5, split="validation")
    item = dataset[0]
    assert item["volume"].shape == (1, 57, 20, 20)
    assert item["volume"].dtype == torch.float32
    assert item["sample_type"].shape == ()
    assert item["sample_type"].dtype == torch.int64
    assert item["slope_xy"].shape == (2,)
    assert item["slope_xy"].dtype == torch.float32
    assert item["metadata"]["source_file"] == "synthetic.root"
    assert item["metadata"]["presence"].dtype == torch.float32
    assert item["metadata"]["hdf5_index"].dtype == torch.int64
    assert item["metadata"]["z_scale"].dtype == torch.float32


def test_validation_uses_exact_central_crop(tmp_path: Path) -> None:
    path = tmp_path / "central.h5"
    y = np.arange(32, dtype=np.int16)[None, :, None]
    x = np.arange(32, dtype=np.int16)[None, None, :]
    z = np.arange(57, dtype=np.int16)[:, None, None]
    one = z * 1000 + y * 32 + x + 2
    regions = np.repeat(one[None], 9, axis=0)
    background = np.ones((9, 57), dtype=np.float32)
    write_hdf5(path, regions=regions, background_mu=background)

    item = HDF5CNN21DDataset(path, split="validation")[0]
    expected_raw = regions[5, :, 6:26, 6:26]
    expected = normalize_per_slice(expected_raw, background[5, 0])
    np.testing.assert_array_equal(item["volume"][0].numpy(), expected)
    np.testing.assert_array_equal(item["metadata"]["crop_offset_xy"].numpy(), [0, 0])


def test_jitter_limits_and_positive_reduction(synthetic_hdf5: Path) -> None:
    config = CNNDataConfig(
        max_jitter=5,
        positive_max_jitter=2,
        random_rotations=False,
        reflect_x=False,
        reflect_y=False,
        seed=77,
    )
    dataset = HDF5CNN21DDataset(
        synthetic_hdf5, split="train", training=True, config=config
    )
    positive_offsets = []
    negative_offsets = []
    for epoch in range(100):
        dataset.set_epoch(epoch)
        positive_offsets.append(dataset[0]["metadata"]["crop_offset_xy"].numpy())
        negative_offsets.append(dataset[1]["metadata"]["crop_offset_xy"].numpy())
    positive_offsets = np.asarray(positive_offsets)
    negative_offsets = np.asarray(negative_offsets)
    assert np.max(np.abs(positive_offsets)) <= 2
    assert np.max(np.abs(negative_offsets)) <= 5
    assert np.max(np.abs(positive_offsets)) == 2
    assert np.max(np.abs(negative_offsets)) == 5


def test_crop_is_plain_slicing_without_wrapping() -> None:
    volume = np.zeros((1, 32, 32), dtype=np.int16)
    volume[0, 10, 31] = 9
    left_crop = extract_xy_crop(volume, crop_size=20, offset_xy=(-5, 0))
    assert left_crop.shape == (1, 20, 20)
    assert left_crop.sum() == 0

    volume.fill(0)
    volume[0, 10, 0] = 11
    right_crop = extract_xy_crop(volume, crop_size=20, offset_xy=(5, 0))
    assert right_crop.sum() == 0
    with pytest.raises(ValueError, match="outside source shape"):
        extract_xy_crop(volume, crop_size=20, offset_xy=(7, 0))


def test_scalar_background_normalization() -> None:
    crop = np.asarray([np.full((2, 3), 5), np.full((2, 3), 13)], dtype=np.int16)
    mu = np.float32(4.0)
    normalized = normalize_per_slice(crop, mu)
    expected = np.asarray(
        [
            np.full((2, 3), 1.0 / np.sqrt(4.0 + 1.0e-6)),
            np.full((2, 3), 9.0 / np.sqrt(4.0 + 1.0e-6)),
        ],
        dtype=np.float32,
    )
    np.testing.assert_allclose(normalized, expected, rtol=1e-6)


def test_local_mpv_normalizations_use_per_layer_reference() -> None:
    volume = np.asarray(
        [
            [[2, 2, 2, 1], [2, 2, 3, 12]],
            [[5, 5, 4, 5], [5, 6, 5, 15]],
        ],
        dtype=np.int16,
    )
    mpv, sigma = local_mpv_sigma_per_z(volume)
    np.testing.assert_array_equal(mpv, [2.0, 5.0])
    assert np.all(sigma > 0.0)
    residual = normalize_per_slice(
        volume, np.float32(1.0), mode="mpv_z_residual", mpv_z=mpv
    )
    np.testing.assert_array_equal(residual[:, 0, :2], 0.0)
    normalized = normalize_per_slice(
        volume, np.float32(1.0), mode="mpv_z_sigma", mpv_z=mpv, sigma_z=sigma
    )
    np.testing.assert_allclose(normalized[:, 0, :2], 0.0)


@pytest.mark.parametrize(
    ("turns", "expected_slope", "marker_yx"),
    [
        (0, (2.0, 3.0), (2, 4)),
        (1, (-3.0, 2.0), (4, 2)),
        (2, (-2.0, -3.0), (2, 0)),
        (3, (3.0, -2.0), (0, 2)),
    ],
)
def test_physical_ccw_rotations_transform_image_and_slope(
    turns: int, expected_slope: tuple[float, float], marker_yx: tuple[int, int]
) -> None:
    image = np.zeros((1, 5, 5), dtype=np.float32)
    image[0, 2, 4] = 1.0  # Synthetic +x displacement in ROOT xy.
    transformed, slope = transform_xy(
        image, np.asarray([2.0, 3.0]), quarter_turns_ccw=turns
    )
    assert transformed[(0, *marker_yx)] == 1.0
    np.testing.assert_array_equal(slope, expected_slope)


@pytest.mark.parametrize(
    ("reflect_x", "reflect_y", "expected_slope", "marker_yx"),
    [
        (True, False, (-2.0, 3.0), (3, 0)),
        (False, True, (2.0, -3.0), (1, 4)),
        (True, True, (-2.0, -3.0), (1, 0)),
    ],
)
def test_reflections_transform_image_and_slope(
    reflect_x: bool,
    reflect_y: bool,
    expected_slope: tuple[float, float],
    marker_yx: tuple[int, int],
) -> None:
    image = np.zeros((1, 5, 5), dtype=np.float32)
    image[0, 3, 4] = 1.0
    transformed, slope = transform_xy(
        image,
        np.asarray([2.0, 3.0]),
        reflect_x=reflect_x,
        reflect_y=reflect_y,
    )
    assert transformed[(0, *marker_yx)] == 1.0
    np.testing.assert_array_equal(slope, expected_slope)


@pytest.mark.parametrize("split", ["validation", "test"])
def test_validation_and_test_are_deterministic(synthetic_hdf5: Path, split: str) -> None:
    dataset = HDF5CNN21DDataset(
        synthetic_hdf5,
        split=split,  # type: ignore[arg-type]
        config=CNNDataConfig(max_jitter=5, seed=1),
    )
    first = dataset[0]
    dataset.set_epoch(99)
    second = dataset[0]
    torch.testing.assert_close(first["volume"], second["volume"], rtol=0, atol=0)
    torch.testing.assert_close(first["slope_xy"], second["slope_xy"], rtol=0, atol=0)
    np.testing.assert_array_equal(second["metadata"]["crop_offset_xy"].numpy(), [0, 0])
    assert int(second["metadata"]["rotation_quarter_turns_ccw"]) == 0
    assert not bool(second["metadata"]["reflected_x"])
    assert not bool(second["metadata"]["reflected_y"])


def test_fixed_seed_reproduces_training_augmentation(synthetic_hdf5: Path) -> None:
    config = CNNDataConfig(seed=9876)
    first = HDF5CNN21DDataset(synthetic_hdf5, split="train", config=config)
    second = HDF5CNN21DDataset(synthetic_hdf5, split="train", config=config)
    for epoch in range(5):
        first.set_epoch(epoch)
        second.set_epoch(epoch)
        one = first[2]
        two = second[2]
        torch.testing.assert_close(one["volume"], two["volume"], rtol=0, atol=0)
        torch.testing.assert_close(one["slope_xy"], two["slope_xy"], rtol=0, atol=0)
        for key in (
            "crop_offset_xy",
            "rotation_quarter_turns_ccw",
            "reflected_x",
            "reflected_y",
        ):
            torch.testing.assert_close(one["metadata"][key], two["metadata"][key])


def test_poisson_background_augmentation_is_train_only_and_updates_mu(
    tmp_path: Path,
) -> None:
    synthetic_hdf5 = tmp_path / "poisson_augmentation.h5"
    regions = np.full((9, *SOURCE_SHAPE), 5, dtype=np.int16)
    regions[:, :, 9:11, 12:14] = 0
    write_hdf5(synthetic_hdf5, regions=regions)
    config = CNNDataConfig(
        max_jitter=0,
        random_rotations=False,
        reflect_x=False,
        reflect_y=False,
        poisson_background_augmentation_probability=1.0,
        poisson_background_mu_ranges=((8.0, 8.0),),
        seed=91,
    )
    train = HDF5CNN21DDataset(
        synthetic_hdf5, split="train", training=True, config=config
    )
    item = train[0]
    assert float(item["metadata"]["background_mu_original"]) == 4.0
    assert float(item["metadata"]["background_mu"]) == 8.0
    assert float(item["metadata"]["poisson_background_added_mu"]) == 4.0
    assert bool(item["metadata"]["poisson_background_augmented"])
    reconstructed = item["volume"][0].numpy() * np.sqrt(8.0 + 1.0e-6) + 8.0
    original = np.asarray(train._file()["regions_raw"][0, :, 6:26, 6:26])
    assert np.all(reconstructed >= original - 1.0e-4)
    np.testing.assert_allclose(reconstructed[original == 0], 0.0, rtol=0, atol=1.0e-4)
    assert np.any(reconstructed[original > 0] > original[original > 0] + 0.5)

    validation = HDF5CNN21DDataset(
        synthetic_hdf5, split="validation", training=False, config=config
    )[0]
    assert float(validation["metadata"]["background_mu"]) == 4.0
    assert float(validation["metadata"]["poisson_background_added_mu"]) == 0.0
    assert not bool(validation["metadata"]["poisson_background_augmented"])


def test_poisson_background_evaluation_target_is_fixed_and_preserves_zeros(
    tmp_path: Path,
) -> None:
    path = tmp_path / "poisson_evaluation.h5"
    regions = np.full((9, *SOURCE_SHAPE), 5, dtype=np.int16)
    regions[:, :, 9:11, 12:14] = 0
    write_hdf5(path, regions=regions)
    config = CNNDataConfig(
        max_jitter=0,
        random_rotations=False,
        reflect_x=False,
        reflect_y=False,
        poisson_background_evaluation_target_mu=7.5,
        seed=92,
    )
    dataset = HDF5CNN21DDataset(path, split="validation", training=False, config=config)
    first = dataset[0]
    second = dataset[0]
    torch.testing.assert_close(first["volume"], second["volume"], rtol=0, atol=0)
    assert float(first["metadata"]["background_mu_original"]) == 4.0
    assert float(first["metadata"]["background_mu"]) == 7.5
    assert float(first["metadata"]["poisson_background_added_mu"]) == 3.5
    reconstructed = first["volume"][0].numpy() * np.sqrt(7.5 + 1.0e-6) + 7.5
    original = np.asarray(dataset._file()["regions_raw"][5, :, 6:26, 6:26])
    np.testing.assert_allclose(reconstructed[original == 0], 0.0, rtol=0, atol=1.0e-4)


def test_layer_drop_augmentation_is_train_only_and_uses_raw_zero_value(
    synthetic_hdf5: Path,
) -> None:
    config = CNNDataConfig(
        max_jitter=0,
        random_rotations=False,
        reflect_x=False,
        reflect_y=False,
        layer_drop_two_consecutive_probability=1.0,
        seed=31,
    )
    train = HDF5CNN21DDataset(
        synthetic_hdf5, split="train", training=True, config=config
    )
    item = train[0]
    assert int(item["metadata"]["layer_drop_count"]) == 2
    first = int(item["metadata"]["layer_drop_first"])
    assert 0 <= first < 56
    expected_missing = -4.0 / np.sqrt(4.0 + 1.0e-6)
    np.testing.assert_allclose(item["volume"][0, first:first + 2], expected_missing)

    validation = HDF5CNN21DDataset(
        synthetic_hdf5, split="validation", training=False, config=config
    )[0]
    assert int(validation["metadata"]["layer_drop_count"]) == 0
    assert int(validation["metadata"]["layer_drop_first"]) == -1


def test_nonconstant_background_over_z_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "invalid_background.h5"
    background = np.full((9, 57), 4.0, dtype=np.float32)
    background[3, 20] = 5.0
    write_hdf5(path, background_mu=background)
    with pytest.raises(ValueError, match="physically scalar per event"):
        HDF5CNN21DDataset(path, split="train")


def test_visibility_audit_reports_risk_and_reduced_jitter(tmp_path: Path) -> None:
    path = tmp_path / "visibility.h5"
    regions = np.ones((9, *SOURCE_SHAPE), dtype=np.int16)
    presence = np.asarray([1, 0, 1, 0, 1, 1, 0, 1, 0], dtype=np.uint8)
    # Give every positive a compact, centrally visible excess.
    regions[presence == 1, :, 12:15, 12:15] = 30
    # A positive excess on the lower/left edge of the central crop disappears
    # for the (+5,+5) translation.
    regions[0] = 1
    regions[0, :, 6:9, 6:9] = 30
    background = np.ones((9, 57), dtype=np.float32)
    write_hdf5(path, regions=regions, background_mu=background, presence=presence)

    risky = HDF5CNN21DDataset(
        path,
        split="train",
        config=CNNDataConfig(max_jitter=5, min_positive_visible_fraction=0.7),
    ).audit_positive_visibility()
    safe = HDF5CNN21DDataset(
        path,
        split="train",
        config=CNNDataConfig(
            max_jitter=5,
            positive_max_jitter=0,
            min_positive_visible_fraction=0.7,
        ),
    ).audit_positive_visibility()
    assert risky["at_risk_examples"] >= 1
    assert safe["at_risk_examples"] == 0
    assert risky["sample_type_labels_changed"] == 0


def test_dataloader_batches_and_worker_local_hdf5_handles(synthetic_hdf5: Path) -> None:
    loaders = create_cnn_dataloaders(
        synthetic_hdf5,
        data_config=CNNDataConfig(seed=44),
        loader_config=CNNLoaderConfig(
            batch_size=2,
            num_workers=2,
            shuffle_train=False,
            seed=44,
        ),
    )
    iterator = iter(loaders["train"])
    batch = next(iterator)
    second_batch = next(iterator)
    assert batch["volume"].shape == (2, 1, 57, 20, 20)
    assert batch["sample_type"].shape == (2,)
    assert batch["slope_xy"].shape == (2, 2)
    worker_ids = set(batch["metadata"]["worker_id"].tolist()) | set(
        second_batch["metadata"]["worker_id"].tolist()
    )
    worker_pids = set(batch["metadata"]["worker_pid"].tolist()) | set(
        second_batch["metadata"]["worker_pid"].tolist()
    )
    assert worker_ids == {0, 1}
    assert len(worker_pids) == 2
    assert 0 not in worker_pids


def test_fixed_loader_seed_reproduces_shuffle_and_augmentation(
    synthetic_hdf5: Path,
) -> None:
    kwargs = dict(
        path=synthetic_hdf5,
        data_config=CNNDataConfig(seed=123),
        loader_config=CNNLoaderConfig(batch_size=3, num_workers=0, seed=456),
    )
    batch_a = next(iter(create_cnn_dataloaders(**kwargs)["train"]))
    batch_b = next(iter(create_cnn_dataloaders(**kwargs)["train"]))
    torch.testing.assert_close(batch_a["volume"], batch_b["volume"], rtol=0, atol=0)
    torch.testing.assert_close(
        batch_a["metadata"]["hdf5_index"], batch_b["metadata"]["hdf5_index"]
    )
