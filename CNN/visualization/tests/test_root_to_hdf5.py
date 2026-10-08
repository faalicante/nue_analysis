from __future__ import annotations

import csv
import json
from pathlib import Path

import h5py
import numpy as np
import pytest
import uproot

import data_preparation.root_to_hdf5 as root_to_hdf5
SHAPE = (57, 32, 32)


def write_root(
    path: Path,
    *,
    counts: np.ndarray,
    background_mu: np.ndarray,
    presence: np.ndarray,
    slope_x: np.ndarray,
    slope_y: np.ndarray,
    signal_event_id: np.ndarray,
    cell_id: np.ndarray,
    tag_cell_id: np.ndarray,
    sample_type: np.ndarray | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    branches = {
        "counts": counts,
        "background_mu": background_mu,
        "presence": presence,
        "slope_x": slope_x,
        "slope_y": slope_y,
        "signal_event_id": signal_event_id,
        "cell_id": cell_id,
        "tag_cell_id": tag_cell_id,
    }
    if sample_type is not None:
        branches["sample_type"] = sample_type
    with uproot.recreate(path) as root_file:
        # mktree with arrays writes an actual TTree.  The assignment shorthand in
        # recent uproot versions writes an RNTuple instead, which is not the source
        # format exercised by this converter.
        root_file.mktree("samples", branches)


def base_columns(count: int) -> dict[str, np.ndarray]:
    return {
        "background_mu": np.full(count, 3.5, dtype=np.float32),
        "presence": np.ones(count, dtype=np.int32),
        "slope_x": np.full(count, 10.0, dtype=np.float32),
        "slope_y": np.full(count, -2.0, dtype=np.float32),
        "signal_event_id": np.arange(count, dtype=np.int32),
        "cell_id": np.full(count, -1, dtype=np.int32),
        "tag_cell_id": np.full(count, -1, dtype=np.int32),
    }


def run_convert(source: Path, output: Path, *extra: str) -> int:
    return root_to_hdf5.main(
        [
            "convert",
            str(source),
            "--output",
            str(output),
            "--chunk-size",
            "2",
            *extra,
        ]
    )


def test_end_to_end_preserves_raw_counts_and_groups(tmp_path: Path) -> None:
    source = tmp_path / "samples.root"
    output = tmp_path / "dataset.h5"
    count = 6
    counts = np.arange(count * np.prod(SHAPE), dtype=np.int32).reshape(count, *SHAPE) % 17
    write_root(
        source,
        counts=counts,
        background_mu=np.asarray([3.0, 3.0, 4.0, 4.0, 5.0, 5.0], dtype=np.float32),
        # Hard-negative presence=1 marks a shower; sample_type distinguishes it from signal.
        presence=np.asarray([1, 1, 0, 1, 0, 1], dtype=np.int32),
        sample_type=np.asarray([2, 2, 0, 1, 0, 1], dtype=np.int32),
        slope_x=np.asarray([10.0, 12.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32),
        slope_y=np.asarray([-2.0, 3.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32),
        signal_event_id=np.asarray([7, 7, -1, -1, -1, -1], dtype=np.int32),
        cell_id=np.asarray([-1, -1, 100, 100, 101, 101], dtype=np.int32),
        tag_cell_id=np.asarray([-1, -1, 1, 2, 1, 2], dtype=np.int32),
    )

    assert run_convert(source, output, "--qa-per-class", "1") == 0

    with h5py.File(output, "r") as hdf5:
        np.testing.assert_array_equal(hdf5["regions_raw"][:], counts)
        assert hdf5["regions_raw"].shape == (count, *SHAPE)
        assert hdf5["regions_raw"].chunks == (2, *SHAPE)
        assert hdf5["regions_raw"].compression == "gzip"
        assert hdf5["background_mu"].shape == (count, 57)
        np.testing.assert_array_equal(hdf5["presence"][:], [1, 1, 0, 1, 0, 1])
        np.testing.assert_array_equal(hdf5["sample_type"][:], [2, 2, 0, 1, 0, 1])
        np.testing.assert_array_equal(hdf5["negative_type"][:], [-1, -1, 0, 1, 0, 1])
        np.testing.assert_allclose(hdf5["slope_xy"][0], [0.27, -0.054], rtol=1e-6)
        split = hdf5["split"][:]
        assert split[0] == split[1]  # signal_event_id=7
        assert split[2] == split[3]  # tile_id=100
        assert split[4] == split[5]  # tile_id=101
        assert set(hdf5["split_indices"].keys()) == {"train", "validation", "test"}

    manifest = output.with_suffix(".manifest.csv")
    with manifest.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == count
    assert rows[0]["group_kind"] == "signal_event_id"
    assert rows[2]["group_kind"] == "tile_id"
    assert list((tmp_path / "dataset_qa").glob("signal_*.png"))
    assert list((tmp_path / "dataset_qa").glob("hard_*.png"))
    assert list((tmp_path / "dataset_qa").glob("poisson_*.png"))


def test_common_background_is_stored_once(tmp_path: Path) -> None:
    source = tmp_path / "common.root"
    output = tmp_path / "common.h5"
    columns = base_columns(2)
    write_root(source, counts=np.ones((2, *SHAPE), dtype=np.int16), **columns)

    assert run_convert(source, output, "--qa-per-class", "0") == 0
    with h5py.File(output, "r") as hdf5:
        assert hdf5["background_mu"].shape == (57,)
        assert bool(hdf5["background_mu"].attrs["common_to_all_events"])


def test_unusable_signal_crop_is_excluded_but_background_is_preserved(tmp_path: Path) -> None:
    source = tmp_path / "support.root"
    output = tmp_path / "support.h5"
    counts = np.ones((3, *SHAPE), dtype=np.int16)
    counts[1, :, :, 19:] = 0  # Only 19 occupied x columns in every slice.
    counts[2] = 0  # Background classes are not subject to the signal-only filter.
    columns = base_columns(3)
    columns.update(
        presence=np.asarray([1, 1, 0], dtype=np.int32),
        sample_type=np.asarray([2, 2, 0], dtype=np.int32),
        signal_event_id=np.asarray([10, 11, -1], dtype=np.int32),
        slope_x=np.asarray([10.0, 10.0, 0.0], dtype=np.float32),
        slope_y=np.asarray([1.0, 1.0, 0.0], dtype=np.float32),
        cell_id=np.asarray([-1, -1, 42], dtype=np.int32),
        tag_cell_id=np.asarray([-1, -1, 1], dtype=np.int32),
    )
    write_root(source, counts=counts, **columns)

    assert run_convert(source, output, "--qa-per-class", "0") == 0
    with h5py.File(output, "r") as hdf5:
        assert hdf5["regions_raw"].shape == (2, *SHAPE)
        np.testing.assert_array_equal(hdf5["source_entry"][:], [0, 2])
        np.testing.assert_array_equal(hdf5["signal_event_id"][:], [10, -1])
        assert int(hdf5.attrs["minimum_signal_support_xy"]) == 20
        assert int(hdf5.attrs["excluded_unusable_signal_entries"]) == 1

    statistics = json.loads(output.with_suffix(".statistics.json").read_text())
    filtering = statistics["signal_support_filter"]
    assert filtering["excluded_signal_entries"] == 1
    assert filtering["excluded_signal_event_ids"] == [11]
    assert filtering["excluded"][0]["minimum_nonzero_width"] == 19


def test_known_bad_background_tile_excludes_hard_and_poisson(tmp_path: Path) -> None:
    source = tmp_path / "hard_tile.root"
    output = tmp_path / "hard_tile.h5"
    counts = np.ones((4, *SHAPE), dtype=np.int16)
    columns = base_columns(4)
    columns.update(
        presence=np.asarray([1, 1, 0, 0], dtype=np.int32),
        sample_type=np.asarray([1, 1, 0, 0], dtype=np.int32),
        signal_event_id=np.full(4, -1, dtype=np.int32),
        slope_x=np.zeros(4, dtype=np.float32),
        slope_y=np.zeros(4, dtype=np.float32),
        cell_id=np.asarray([58, 59, 58, 59], dtype=np.int32),
        tag_cell_id=np.asarray([2, 3, 4, 5], dtype=np.int32),
    )
    write_root(source, counts=counts, **columns)

    assert run_convert(source, output, "--qa-per-class", "0") == 0
    with h5py.File(output, "r") as hdf5:
        np.testing.assert_array_equal(hdf5["source_entry"][:], [1, 3])
        np.testing.assert_array_equal(hdf5["sample_type"][:], [1, 0])
        np.testing.assert_array_equal(hdf5["tile_id"][:], [59, 59])
        assert int(hdf5.attrs["excluded_background_entries"]) == 2

    statistics = json.loads(output.with_suffix(".statistics.json").read_text())
    filtering = statistics["background_tile_filter"]
    assert filtering["excluded_background_tile_ids"] == [58]
    assert filtering["excluded_background_entries"] == 2
    assert filtering["excluded_by_class"] == {"poisson": 1, "hard": 1}
    assert filtering["excluded"][0]["crop_id_in_tile"] == 2


def test_dimension_mismatch_stops_without_output(tmp_path: Path, capsys) -> None:
    source = tmp_path / "wrong.root"
    output = tmp_path / "wrong.h5"
    columns = base_columns(1)
    write_root(source, counts=np.zeros((1, 57, 20, 20), dtype=np.int32), **columns)

    assert run_convert(source, output, "--qa-per-class", "0") == 2
    assert not output.exists()
    assert "Dimension mismatch" in capsys.readouterr().err


def test_invalid_values_stop_conversion(tmp_path: Path, capsys) -> None:
    source = tmp_path / "invalid.root"
    output = tmp_path / "invalid.h5"
    columns = base_columns(1)
    columns["background_mu"][0] = 0.0
    write_root(source, counts=np.zeros((1, *SHAPE), dtype=np.int32), **columns)

    assert run_convert(source, output, "--qa-per-class", "0") == 2
    assert not output.exists()
    assert "background_mu <= 0" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("negative_counts", "Negative raw counts"),
        ("nan_slope", "NaN or infinite slope"),
        ("infinite_background", "NaN or infinite background_mu"),
        ("invalid_label", "Invalid labels"),
    ],
)
def test_other_invalid_physical_values_stop(
    tmp_path: Path, capsys, case: str, message: str
) -> None:
    source = tmp_path / f"{case}.root"
    output = tmp_path / f"{case}.h5"
    counts = np.zeros((1, *SHAPE), dtype=np.int32)
    columns = base_columns(1)
    if case == "negative_counts":
        counts[0, 0, 0, 0] = -1
    elif case == "nan_slope":
        columns["slope_x"][0] = np.nan
    elif case == "infinite_background":
        columns["background_mu"][0] = np.inf
    elif case == "invalid_label":
        columns["presence"][0] = 3
    write_root(source, counts=counts, **columns)

    assert run_convert(source, output, "--qa-per-class", "0") == 2
    assert not output.exists()
    assert message in capsys.readouterr().err


def test_existing_splits_are_preserved_and_leakage_is_rejected(tmp_path: Path, capsys) -> None:
    train_source = tmp_path / "train" / "background.root"
    validation_source = tmp_path / "validation" / "background.root"
    preserved_output = tmp_path / "preserved.h5"
    leak_output = tmp_path / "leak.h5"
    columns = base_columns(1)
    columns.update(
        presence=np.zeros(1, dtype=np.int32),
        slope_x=np.zeros(1, dtype=np.float32),
        slope_y=np.zeros(1, dtype=np.float32),
        signal_event_id=np.full(1, -1, dtype=np.int32),
        cell_id=np.asarray([42], dtype=np.int32),
        tag_cell_id=np.asarray([1], dtype=np.int32),
    )
    counts = np.zeros((1, *SHAPE), dtype=np.int32)
    write_root(train_source, counts=counts, **columns)
    columns["cell_id"][0] = 43
    write_root(validation_source, counts=counts, **columns)

    assert (
        root_to_hdf5.main(
            [
                "convert",
                str(train_source),
                str(validation_source),
                "--output",
                str(preserved_output),
                "--qa-per-class",
                "0",
            ]
        )
        == 0
    )
    with h5py.File(preserved_output, "r") as hdf5:
        np.testing.assert_array_equal(hdf5["split"][:], [0, 1])

    # Re-create the validation source with the same tile as train: preserving the
    # existing directory split would now leak tile_id=42 and must be rejected.
    columns["cell_id"][0] = 42
    write_root(validation_source, counts=counts, **columns)

    assert (
        root_to_hdf5.main(
            [
                "convert",
                str(train_source),
                str(validation_source),
                "--output",
                str(leak_output),
                "--qa-per-class",
                "0",
            ]
        )
        == 2
    )
    assert not leak_output.exists()
    assert "Existing split leakage" in capsys.readouterr().err
