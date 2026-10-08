from __future__ import annotations

import h5py
import numpy as np
import uproot

from scanning.produce_background_scan_shards import selected_cells
from scanning.root_background_xypseg_to_scan_hdf5 import (
    BackgroundCell,
    _assert_physical_cell_edges,
    render_background_root_path,
    validate_background_root,
    write_background_hdf5,
    write_partial_background_hdf5,
)
from scanning.root_xypseg_to_scan_hdf5 import extract_volume


def test_background_cell_geometry_and_path() -> None:
    cell = BackgroundCell(1, 1)
    assert cell.cell_id == 0
    assert cell.x_bounds_um == (200_000.0, 210_000.0)
    assert cell.y_bounds_um == (4_500.0, 14_500.0)
    assert cell.center_um == (205_000.0, 9_500.0)
    path = render_background_root_path(
        "/cell_{x}0_{y}0/b000021.0.{x}.{y}.{cell_id}.root", cell
    )
    assert path == "/cell_10_10/b000021.0.1.1.0.root"


def test_background_cell_custom_grid_origin() -> None:
    first = BackgroundCell(1, 1, grid_xmin_um=5_729.0, grid_ymin_um=196_206.0)
    second = BackgroundCell(2, 3, grid_xmin_um=5_729.0, grid_ymin_um=196_206.0)
    assert first.x_bounds_um == (5_729.0, 15_729.0)
    assert first.y_bounds_um == (196_206.0, 206_206.0)
    assert first.center_um == (10_729.0, 201_206.0)
    assert second.x_bounds_um == (15_729.0, 25_729.0)
    assert second.y_bounds_um == (216_206.0, 226_206.0)


def test_cell_id_order_and_known_bad_exclusion() -> None:
    cells = selected_cells(1, 18, 1, 18, include_known_bad=False)
    assert len(cells) == 323
    assert [cell.cell_id for cell in cells[:3]] == [0, 1, 2]
    assert 58 not in {cell.cell_id for cell in cells}
    included = selected_cells(1, 18, 1, 18, include_known_bad=True)
    assert len(included) == 324
    bad = next(cell for cell in included if cell.cell_id == 58)
    assert (bad.x, bad.y) == (5, 4)


def test_realistic_599_bin_axis_accepts_sub_bin_boundary_offset() -> None:
    cell = BackgroundCell(9, 9)
    x_edges = np.linspace(270_000.3125, 299_999.65625, 600)
    y_edges = np.linspace(74_500.1484375, 104_499.9140625, 600)
    values = np.zeros((599, 599, 60), dtype=np.uint8)
    _, local_x_edges, local_y_edges, _ = extract_volume(
        values,
        x_edges,
        y_edges,
        center_x_um=cell.center_um[0],
        center_y_um=cell.center_um[1],
        size=200,
        plates=57,
    )
    _assert_physical_cell_edges(cell, local_x_edges, local_y_edges)
    assert abs(local_x_edges[0] - cell.x_bounds_um[0]) < 50.1
    assert abs(local_x_edges[-1] - cell.x_bounds_um[1]) < 50.1
    assert abs(local_y_edges[0] - cell.y_bounds_um[0]) < 50.1
    assert abs(local_y_edges[-1] - cell.y_bounds_um[1]) < 50.1


def test_data_brick_122_complete_map_accepts_shifted_axis_origin() -> None:
    """A full 200-bin map is valid when its centre and span agree within one bin."""
    from scanning.produce_data_scan_shards import DataCell

    cell = DataCell(3, 3)
    # Exact Y bounds reported by b000122 cell_30_30.  The lower endpoint is
    # 68.5 um from nominal, but the centre and 200-bin span are each correct to
    # less than the 50.228 um source-bin width.
    y_edges = np.linspace(24_931.50477636256, 34_977.135700533174, 201)
    x_edges = np.linspace(24_954.0, 34_999.6, 201)
    _assert_physical_cell_edges(cell, x_edges, y_edges)


def test_shifted_map_still_rejects_the_adjacent_cell() -> None:
    from scanning.produce_data_scan_shards import DataCell

    cell = DataCell(4, 3)
    edges = np.linspace(24_931.50477636256, 34_977.135700533174, 201)
    with np.testing.assert_raises_regex(ValueError, "center error"):
        _assert_physical_cell_edges(cell, edges, edges)


def test_write_background_hdf5_without_truth(tmp_path) -> None:
    cell = BackgroundCell(1, 1)
    x_edges = np.arange(199_000.0, 211_000.0 + 50.0, 50.0)
    y_edges = np.arange(3_500.0, 15_500.0 + 50.0, 50.0)
    z_edges = np.arange(60, dtype=np.float64)
    values = np.zeros((240, 240, 59), dtype=np.float32)
    values[20:220, 20:220, 1:58] = 4.0
    xyseg = np.full((240, 240), 197.0, dtype=np.float32)
    root_path = tmp_path / "cell.root"
    with uproot.recreate(root_path) as root_file:
        root_file["XYPseg"] = values, x_edges, y_edges, z_edges
        root_file["XYseg"] = xyseg, x_edges, y_edges

    validate_background_root(str(root_path), cell)
    output = tmp_path / "background.h5"
    report = write_background_hdf5([(cell, str(root_path))], output)
    assert report["size_bytes"] > 0
    with h5py.File(output, "r") as hdf5:
        assert "truth" not in hdf5
        assert hdf5["volumes_raw"].shape == (1, 57, 200, 200)
        assert int(hdf5["volumes_raw"][:].sum()) == 57 * 200 * 200 * 4
        assert int(hdf5["event_id"][0]) == 0
        assert int(hdf5["cell_x"][0]) == 1
        assert int(hdf5["cell_y"][0]) == 1
        np.testing.assert_allclose(hdf5["x_edges_um"][0, [0, -1]], [200_000, 210_000])
        np.testing.assert_allclose(hdf5["y_edges_um"][0, [0, -1]], [4_500, 14_500])


def test_validate_background_root_rejects_missing_background_estimate(tmp_path) -> None:
    cell = BackgroundCell(1, 1)
    x_edges = np.arange(199_000.0, 211_000.0 + 50.0, 50.0)
    y_edges = np.arange(3_500.0, 15_500.0 + 50.0, 50.0)
    z_edges = np.arange(60, dtype=np.float64)
    values = np.zeros((240, 240, 59), dtype=np.float32)
    values[20:220, 20:220, 1:58] = 4.0
    root_path = tmp_path / "empty_xyseg.root"
    with uproot.recreate(root_path) as root_file:
        root_file["XYPseg"] = values, x_edges, y_edges, z_edges
        root_file["XYseg"] = np.zeros((240, 240), dtype=np.float32), x_edges, y_edges

    with np.testing.assert_raises_regex(ValueError, "XYseg has no bins with content > 50"):
        validate_background_root(str(root_path), cell)


def test_write_partial_background_hdf5_without_padding(tmp_path) -> None:
    cell = BackgroundCell(1, 1)
    x_edges = np.arange(200_000.0, 210_000.0 + 50.0, 50.0)
    y_edges = np.arange(4_500.0, 9_500.0 + 50.0, 50.0)
    z_edges = np.arange(60, dtype=np.float64)
    values = np.zeros((200, 100, 59), dtype=np.float32)
    values[:, :, 1:58] = 4.0
    xyseg = np.full((200, 100), 197.0, dtype=np.float32)
    root_path = tmp_path / "partial.root"
    with uproot.recreate(root_path) as root_file:
        root_file["XYPseg"] = values, x_edges, y_edges, z_edges
        root_file["XYseg"] = xyseg, x_edges, y_edges

    output = tmp_path / "partial.h5"
    report = write_partial_background_hdf5(cell, str(root_path), output)
    assert report["partial_cell"] is True
    with h5py.File(output, "r") as hdf5:
        assert hdf5["volumes_raw"].shape == (1, 57, 100, 200)
        assert hdf5.attrs["partial_cell"]
        np.testing.assert_allclose(hdf5["y_edges_um"][0, [0, -1]], [4_500, 9_500])
