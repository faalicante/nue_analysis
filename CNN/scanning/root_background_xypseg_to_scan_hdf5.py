#!/usr/bin/env python3
"""Extract physical 57x200x200 background cells from XYPseg ROOT files."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import h5py
import numpy as np
import uproot

from scanning.root_xypseg_to_scan_hdf5 import (
    CropSupportError,
    DEFAULT_PLATES,
    DEFAULT_VOLUME_SIZE,
    DEFAULT_XY_BIN_SIZE_UM,
    extract_available_volume,
    extract_volume,
    poisson_background_per_layer,
)


CELL_SIZE_UM = 10_000.0
DEFAULT_CELL_XMIN_UM = 200_000.0
DEFAULT_CELL_YMIN_UM = 4_500.0


@dataclass(frozen=True)
class BackgroundCell:
    """One physical brick cell, using the one-based x/y convention of crop_bw.C."""

    x: int
    y: int
    grid_xmin_um: float = DEFAULT_CELL_XMIN_UM
    grid_ymin_um: float = DEFAULT_CELL_YMIN_UM

    def __post_init__(self) -> None:
        if not 1 <= self.x <= 18 or not 1 <= self.y <= 18:
            raise ValueError(f"cell coordinates must be in [1, 18], got ({self.x}, {self.y})")

    @property
    def cell_id(self) -> int:
        return (self.y - 1) * 18 + (self.x - 1)

    @property
    def x_bounds_um(self) -> tuple[float, float]:
        low = self.grid_xmin_um + (self.x - 1) * CELL_SIZE_UM
        return low, low + CELL_SIZE_UM

    @property
    def y_bounds_um(self) -> tuple[float, float]:
        low = self.grid_ymin_um + (self.y - 1) * CELL_SIZE_UM
        return low, low + CELL_SIZE_UM

    @property
    def center_um(self) -> tuple[float, float]:
        x_low, x_high = self.x_bounds_um
        y_low, y_high = self.y_bounds_um
        return 0.5 * (x_low + x_high), 0.5 * (y_low + y_high)


def render_background_root_path(template: str, cell: BackgroundCell) -> str:
    try:
        return template.format(x=cell.x, y=cell.y, cell_id=cell.cell_id)
    except KeyError as error:
        raise ValueError("root template may only use {x}, {y}, and {cell_id}") from error


def _assert_physical_cell_edges(
    cell: BackgroundCell,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> None:
    expected_x = np.asarray(cell.x_bounds_um)
    expected_y = np.asarray(cell.y_bounds_um)
    observed_x = np.asarray((x_edges[0], x_edges[-1]))
    observed_y = np.asarray((y_edges[0], y_edges[-1]))
    # Real files need not have exactly 50 um bins and their ROOT axis origin can
    # be shifted with respect to the nominal cell grid.  Consequently a complete
    # 200-bin map can differ from either nominal boundary by more than one bin
    # even though its centre and physical span are both correct within one bin.
    # Validate those two invariant quantities instead of each endpoint.
    x_bin_width = float(np.median(np.diff(x_edges)))
    y_bin_width = float(np.median(np.diff(y_edges)))
    if not np.isclose(x_bin_width, DEFAULT_XY_BIN_SIZE_UM, rtol=0.02, atol=0.0):
        raise ValueError(f"unexpected x bin width {x_bin_width} um")
    if not np.isclose(y_bin_width, DEFAULT_XY_BIN_SIZE_UM, rtol=0.02, atol=0.0):
        raise ValueError(f"unexpected y bin width {y_bin_width} um")
    observed_x_center = float(observed_x.mean())
    observed_y_center = float(observed_y.mean())
    expected_x_center = float(expected_x.mean())
    expected_y_center = float(expected_y.mean())
    observed_x_span = float(observed_x[1] - observed_x[0])
    observed_y_span = float(observed_y[1] - observed_y[0])
    expected_x_span = float(expected_x[1] - expected_x[0])
    expected_y_span = float(expected_y[1] - expected_y[0])
    x_matches = (
        abs(observed_x_center - expected_x_center) <= x_bin_width
        and abs(observed_x_span - expected_x_span) <= x_bin_width
    )
    y_matches = (
        abs(observed_y_center - expected_y_center) <= y_bin_width
        and abs(observed_y_span - expected_y_span) <= y_bin_width
    )
    if not x_matches:
        raise ValueError(
            f"cell ({cell.x}, {cell.y}) x edges {observed_x.tolist()} do not match "
            f"physical bounds {expected_x.tolist()}; center error "
            f"{observed_x_center - expected_x_center:.3f} um, span error "
            f"{observed_x_span - expected_x_span:.3f} um, bin width {x_bin_width:.3f} um"
        )
    if not y_matches:
        raise ValueError(
            f"cell ({cell.x}, {cell.y}) y edges {observed_y.tolist()} do not match "
            f"physical bounds {expected_y.tolist()}; center error "
            f"{observed_y_center - expected_y_center:.3f} um, span error "
            f"{observed_y_span - expected_y_span:.3f} um, bin width {y_bin_width:.3f} um"
        )


def validate_background_root(path: str, cell: BackgroundCell) -> None:
    """Validate keys, z support and exact 200x200 physical-cell alignment."""
    with uproot.open(path) as root_file:
        if "XYPseg" not in root_file or "XYseg" not in root_file:
            raise KeyError(f"{path}: expected XYPseg and XYseg")
        histogram = root_file["XYPseg"]
        shape = histogram.values(flow=False).shape
        xyseg = root_file["XYseg"].values(flow=False)
        if len(shape) != 3 or shape[2] < DEFAULT_PLATES + 1:
            raise ValueError(
                f"{path}: XYPseg shape {shape} cannot provide ROOT z bins 2..58"
            )
        x_edges = np.asarray(histogram.axes[0].edges(), dtype=np.float64)
        y_edges = np.asarray(histogram.axes[1].edges(), dtype=np.float64)
        center_x, center_y = cell.center_um
        dummy = np.empty(shape, dtype=np.uint8)
        try:
            _, local_x_edges, local_y_edges, _ = extract_volume(
                dummy,
                x_edges,
                y_edges,
                center_x_um=center_x,
                center_y_um=center_y,
                size=DEFAULT_VOLUME_SIZE,
                plates=DEFAULT_PLATES,
            )
        except ValueError as error:
            raise CropSupportError(str(error)) from error
        _assert_physical_cell_edges(cell, local_x_edges, local_y_edges)
        try:
            poisson_background_per_layer(xyseg)
        except ValueError as error:
            raise ValueError(f"{path}: {error}") from error


def validate_partial_background_root(
    path: str,
    cell: BackgroundCell,
    *,
    minimum_xy_bins: int = 20,
) -> None:
    """Validate a clipped data cell that is large enough for CNN windows."""
    with uproot.open(path) as root_file:
        if "XYPseg" not in root_file or "XYseg" not in root_file:
            raise KeyError(f"{path}: expected XYPseg and XYseg")
        histogram = root_file["XYPseg"]
        values = histogram.values(flow=False)
        xyseg = root_file["XYseg"].values(flow=False)
        x_edges = np.asarray(histogram.axes[0].edges(), dtype=np.float64)
        y_edges = np.asarray(histogram.axes[1].edges(), dtype=np.float64)
    x_bin_width = float(np.median(np.diff(x_edges)))
    y_bin_width = float(np.median(np.diff(y_edges)))
    if not np.isclose(x_bin_width, DEFAULT_XY_BIN_SIZE_UM, rtol=0.02, atol=0.0):
        raise ValueError(f"unexpected x bin width {x_bin_width} um")
    if not np.isclose(y_bin_width, DEFAULT_XY_BIN_SIZE_UM, rtol=0.02, atol=0.0):
        raise ValueError(f"unexpected y bin width {y_bin_width} um")
    extract_available_volume(
        values,
        x_edges,
        y_edges,
        center_x_um=cell.center_um[0],
        center_y_um=cell.center_um[1],
        plates=DEFAULT_PLATES,
        minimum_xy_bins=minimum_xy_bins,
    )
    try:
        poisson_background_per_layer(xyseg)
    except ValueError as error:
        raise ValueError(f"{path}: {error}") from error


def write_background_hdf5(
    cells_and_paths: Sequence[tuple[BackgroundCell, str]],
    output: Path,
) -> dict[str, object]:
    """Write a truth-free scan shard consumable by scanning/scan_cnn21d_volumes.py."""
    output.parent.mkdir(parents=True, exist_ok=True)
    n = len(cells_and_paths)
    if n == 0:
        raise ValueError("cannot write an empty background shard")
    string_dtype = h5py.string_dtype(encoding="utf-8")
    report: dict[str, object] = {"cells": [], "output": str(output)}

    with h5py.File(output, "w") as hdf5:
        volumes = hdf5.create_dataset(
            "volumes_raw",
            shape=(n, DEFAULT_PLATES, DEFAULT_VOLUME_SIZE, DEFAULT_VOLUME_SIZE),
            dtype=np.int32,
            chunks=(1, 1, DEFAULT_VOLUME_SIZE, DEFAULT_VOLUME_SIZE),
            compression="gzip",
            compression_opts=4,
            shuffle=True,
        )
        event_ids = hdf5.create_dataset("event_id", shape=(n,), dtype=np.int64)
        cell_ids = hdf5.create_dataset("cell_id", shape=(n,), dtype=np.int16)
        cell_x = hdf5.create_dataset("cell_x", shape=(n,), dtype=np.int8)
        cell_y = hdf5.create_dataset("cell_y", shape=(n,), dtype=np.int8)
        background_ds = hdf5.create_dataset("background_mu", shape=(n,), dtype=np.float32)
        source_files = hdf5.create_dataset("source_file", shape=(n,), dtype=string_dtype)
        x_edges_ds = hdf5.create_dataset(
            "x_edges_um", shape=(n, DEFAULT_VOLUME_SIZE + 1), dtype=np.float64
        )
        y_edges_ds = hdf5.create_dataset(
            "y_edges_um", shape=(n, DEFAULT_VOLUME_SIZE + 1), dtype=np.float64
        )

        hdf5.attrs["sample_kind"] = "background_full_cell"
        hdf5.attrs["volume_axis_order"] = "event,z,y,x"
        hdf5.attrs["source_histogram"] = "XYPseg"
        hdf5.attrs["source_histogram_axis_order"] = "x,y,plate"
        hdf5.attrs["source_root_z_bins"] = "2..58 inclusive"
        hdf5.attrs["cell_coordinate_convention"] = "one-based x,y in [1,18]"
        hdf5.attrs["cell_id_formula"] = "(y - 1) * 18 + (x - 1)"
        hdf5.attrs["xy_bin_size_nominal_um"] = DEFAULT_XY_BIN_SIZE_UM

        for index, (cell, root_path) in enumerate(cells_and_paths):
            with uproot.open(root_path) as root_file:
                histogram = root_file["XYPseg"]
                values = histogram.values(flow=False)
                x_edges = np.asarray(histogram.axes[0].edges(), dtype=np.float64)
                y_edges = np.asarray(histogram.axes[1].edges(), dtype=np.float64)
                xyseg = root_file["XYseg"].values(flow=False)
            if not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-5):
                raise ValueError(f"non-integer XYPseg contents in {root_path}")
            center_x, center_y = cell.center_um
            volume, local_x_edges, local_y_edges, center_bins = extract_volume(
                values,
                x_edges,
                y_edges,
                center_x_um=center_x,
                center_y_um=center_y,
                size=DEFAULT_VOLUME_SIZE,
                plates=DEFAULT_PLATES,
            )
            _assert_physical_cell_edges(cell, local_x_edges, local_y_edges)
            background_mu = poisson_background_per_layer(xyseg)

            volumes[index] = volume
            event_ids[index] = cell.cell_id
            cell_ids[index] = cell.cell_id
            cell_x[index] = cell.x
            cell_y[index] = cell.y
            background_ds[index] = background_mu
            source_files[index] = root_path
            x_edges_ds[index] = local_x_edges
            y_edges_ds[index] = local_y_edges
            report["cells"].append({  # type: ignore[union-attr]
                "cell_id": cell.cell_id,
                "cell_x": cell.x,
                "cell_y": cell.y,
                "source_file": root_path,
                "source_shape_xyz": list(values.shape),
                "volume_shape_zyx": list(volume.shape),
                "volume_sum": int(volume.sum()),
                "background_mu_per_layer": background_mu,
                "center_root_bins": list(center_bins),
                "x_bounds_um": list(cell.x_bounds_um),
                "y_bounds_um": list(cell.y_bounds_um),
            })

    report["size_bytes"] = output.stat().st_size
    return report


def write_partial_background_hdf5(
    cell: BackgroundCell,
    root_path: str,
    output: Path,
    *,
    minimum_xy_bins: int = 20,
) -> dict[str, object]:
    """Write the real intersection of a nominal cell crop and clipped support."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with uproot.open(root_path) as root_file:
        histogram = root_file["XYPseg"]
        values = histogram.values(flow=False)
        x_edges = np.asarray(histogram.axes[0].edges(), dtype=np.float64)
        y_edges = np.asarray(histogram.axes[1].edges(), dtype=np.float64)
        xyseg = root_file["XYseg"].values(flow=False)
    if not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-5):
        raise ValueError(f"non-integer XYPseg contents in {root_path}")
    volume, local_x_edges, local_y_edges = extract_available_volume(
        values,
        x_edges,
        y_edges,
        center_x_um=cell.center_um[0],
        center_y_um=cell.center_um[1],
        plates=DEFAULT_PLATES,
        minimum_xy_bins=minimum_xy_bins,
    )
    try:
        background_mu = poisson_background_per_layer(xyseg)
    except ValueError as error:
        raise ValueError(f"{root_path}: {error}") from error
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(output, "w") as hdf5:
        hdf5.create_dataset(
            "volumes_raw", data=volume[None], dtype=np.int32,
            chunks=(1, 1, volume.shape[1], volume.shape[2]),
            compression="gzip", compression_opts=4, shuffle=True,
        )
        hdf5.create_dataset("event_id", data=np.asarray([cell.cell_id], dtype=np.int64))
        hdf5.create_dataset("cell_id", data=np.asarray([cell.cell_id], dtype=np.int16))
        hdf5.create_dataset("cell_x", data=np.asarray([cell.x], dtype=np.int8))
        hdf5.create_dataset("cell_y", data=np.asarray([cell.y], dtype=np.int8))
        hdf5.create_dataset("background_mu", data=np.asarray([background_mu], dtype=np.float32))
        hdf5.create_dataset("source_file", data=np.asarray([root_path], dtype=object), dtype=string_dtype)
        hdf5.create_dataset("x_edges_um", data=local_x_edges[None])
        hdf5.create_dataset("y_edges_um", data=local_y_edges[None])
        hdf5.attrs["sample_kind"] = "background_partial_cell"
        hdf5.attrs["volume_axis_order"] = "event,z,y,x"
        hdf5.attrs["source_histogram"] = "XYPseg"
        hdf5.attrs["source_histogram_axis_order"] = "x,y,plate"
        hdf5.attrs["source_root_z_bins"] = "2..58 inclusive"
        hdf5.attrs["cell_coordinate_convention"] = "one-based x,y in [1,18]"
        hdf5.attrs["xy_bin_size_nominal_um"] = DEFAULT_XY_BIN_SIZE_UM
        hdf5.attrs["partial_cell"] = True
        hdf5.attrs["nominal_cell_crop_size_bins"] = DEFAULT_VOLUME_SIZE
        hdf5.attrs["minimum_scan_crop_size_bins"] = minimum_xy_bins
    return {
        "output": str(output),
        "size_bytes": output.stat().st_size,
        "cell_id": cell.cell_id,
        "cell_x": cell.x,
        "cell_y": cell.y,
        "source_file": root_path,
        "source_shape_xyz": list(values.shape),
        "volume_shape_zyx": list(volume.shape),
        "x_bounds_um": [float(local_x_edges[0]), float(local_x_edges[-1])],
        "y_bounds_um": [float(local_y_edges[0]), float(local_y_edges[-1])],
        "background_mu_per_layer": background_mu,
        "partial_cell": True,
    }


def write_report(report: dict[str, object], path: Path) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
