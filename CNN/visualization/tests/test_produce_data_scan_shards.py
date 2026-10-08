from pathlib import Path

import h5py
import numpy as np
import pytest

from scanning.produce_data_scan_shards import (
    DataCell,
    add_data_provenance,
    cell_from_id,
    normalize_rwb,
    read_cell_list,
    render_data_root_path,
)
from scanning.root_background_xypseg_to_scan_hdf5 import BackgroundCell


@pytest.mark.parametrize("value", ["21", "021", "b000021"])
def test_normalize_rwb_equivalent_spellings(value: str) -> None:
    assert normalize_rwb(value) == (21, "021", "b000021")


def test_render_data_root_path() -> None:
    result = render_data_root_path(
        "/eos/data/cell_{x}0_{y}0/{brick}/{brick}.0.0.0.trk.root",
        DataCell(x=5, y=4),
        rwb="021",
        brick="b000021",
    )
    assert result == "/eos/data/cell_50_40/b000021/b000021.0.0.0.trk.root"


def test_cell_from_id() -> None:
    assert cell_from_id(58) == DataCell(x=5, y=4)
    assert cell_from_id(323) == DataCell(x=18, y=18)


def test_data_cell_152_geometry() -> None:
    cell = cell_from_id(152)
    assert cell == DataCell(x=9, y=9)
    assert cell.center_um == (90_000.0, 90_000.0)
    assert cell.x_bounds_um == (85_000.0, 95_000.0)
    assert cell.y_bounds_um == (85_000.0, 95_000.0)


def test_read_cell_list_supported_formats(tmp_path: Path) -> None:
    path = tmp_path / "good_cells.txt"
    path.write_text(
        "# selected cells\n"
        "cell_id x y\n"
        "0 1 1\n"
        "2,1\n"
        "58\n"
        "/eos/production/cell_180_180/b000021.0.0.0.trk.root\n"
        "58  # duplicate\n",
        encoding="utf-8",
    )
    assert read_cell_list(path) == [
        DataCell(1, 1),
        DataCell(2, 1),
        DataCell(5, 4),
        DataCell(18, 18),
    ]


def test_read_cell_list_rejects_inconsistent_columns(tmp_path: Path) -> None:
    path = tmp_path / "bad_cells.txt"
    path.write_text("59 5 4\n", encoding="utf-8")
    with pytest.raises(ValueError, match="disagrees"):
        read_cell_list(path)


def test_add_data_provenance(tmp_path: Path) -> None:
    path = tmp_path / "shard.h5"
    with h5py.File(path, "w") as hdf5:
        hdf5.create_dataset("event_id", data=np.asarray([0, 58, 323]))
        hdf5.create_dataset("cell_id", data=np.asarray([0, 58, 323]))

    add_data_provenance(path, rwb_id=21, rwb="021", brick="b000021")

    with h5py.File(path, "r") as hdf5:
        assert hdf5.attrs["sample_kind"] == "data_full_cell"
        assert hdf5["event_id"][:].tolist() == [21000, 21058, 21323]
        assert hdf5["run_wall_brick"][:].tolist() == [21, 21, 21]
        assert hdf5["brick_name"].asstr()[:].tolist() == ["b000021"] * 3
