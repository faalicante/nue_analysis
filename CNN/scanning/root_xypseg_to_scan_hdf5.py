#!/usr/bin/env python3
"""Extract vertex-centred 57x200x200 scan volumes from signal XYPseg TH3 files."""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import h5py
import numpy as np
import uproot


SLOPE_CORRECTION = 1315.0 / 1350.0
DEFAULT_PLATES = 57
DEFAULT_VOLUME_SIZE = 200
DEFAULT_Z_BIN_SIZE_UM = 1350.0
DEFAULT_XY_BIN_SIZE_UM = 50.0
PLATE_INDEX_OFFSET = 3


class CropSupportError(ValueError):
    """A requested fixed-size volume does not fit inside the source histogram."""


@dataclass(frozen=True)
class SignalTruth:
    event_id: int
    xpos_um: float
    ypos_um: float
    projx_um: float
    projy_um: float
    ltx_mrad_raw: float
    lty_mrad_raw: float
    energy_gev: float
    plate: int

    @property
    def ltx_mrad_corrected(self) -> float:
        return self.ltx_mrad_raw * SLOPE_CORRECTION

    @property
    def lty_mrad_corrected(self) -> float:
        return self.lty_mrad_raw * SLOPE_CORRECTION

    def propagated_center_um(self, plates: int = DEFAULT_PLATES) -> tuple[float, float]:
        # Reproduce data_preparation/loop.py: the plate stored in nue_int_10k.txt is converted
        # to the crop_bw.C convention with ``p0 = plate - 3``.
        crop_plate = self.plate - PLATE_INDEX_OFFSET
        propagated_plates = min(40, plates - crop_plate)
        factor = DEFAULT_Z_BIN_SIZE_UM * 0.5 * propagated_plates / 1000.0
        return (
            self.xpos_um + self.ltx_mrad_corrected * factor,
            self.ypos_um + self.lty_mrad_corrected * factor,
        )


def read_truth_txt(path: Path) -> dict[int, SignalTruth]:
    truth: dict[int, SignalTruth] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, skipinitialspace=True), start=1):
            if not row or all(not field.strip() for field in row):
                continue
            if len(row) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 columns, found {len(row)}")
            event_id = int(row[0])
            if event_id in truth:
                raise ValueError(f"{path}:{line_number}: duplicate event {event_id}")
            truth[event_id] = SignalTruth(
                event_id=event_id,
                xpos_um=float(row[1]), ypos_um=float(row[2]),
                projx_um=float(row[3]), projy_um=float(row[4]),
                ltx_mrad_raw=float(row[5]), lty_mrad_raw=float(row[6]),
                energy_gev=float(row[7]), plate=int(row[8]),
            )
    if not truth:
        raise ValueError(f"no truth rows found in {path}")
    return truth


def event_id_from_signal_filename(path: str) -> int:
    name = path.rsplit("/", 1)[-1]
    match = re.search(r"\.(\d+)\.trk\.root$", name)
    if match is None:
        raise ValueError(f"cannot infer event id from signal filename: {path}")
    # crop_bw.C opens cell/event N as ...trk.root with the final field N+1.
    return int(match.group(1)) - 1


def root_find_fix_bin(edges: np.ndarray, value: float) -> int:
    """Return ROOT-style one-based FindFixBin for an in-range value."""
    if value < edges[0] or value >= edges[-1]:
        raise ValueError(f"coordinate {value} outside [{edges[0]}, {edges[-1]})")
    return int(np.searchsorted(edges, value, side="right"))


def crop_bounds(edges: np.ndarray, center: float, size: int) -> tuple[int, int, int]:
    """Return Python start/stop and one-based ROOT centre bin, matching crop_bw.C."""
    if size <= 0 or size % 2:
        raise ValueError("crop size must be a positive even integer")
    center_bin = root_find_fix_bin(edges, center)
    first_bin = center_bin - size // 2
    start = first_bin - 1
    stop = start + size
    if start < 0 or stop > len(edges) - 1:
        raise CropSupportError(f"crop of size {size} around {center} exceeds histogram support")
    return start, stop, center_bin


def validate_root_volume_support(
    root_path: str,
    truth: SignalTruth,
    size: int = DEFAULT_VOLUME_SIZE,
) -> None:
    """Read only TH3 axes and verify that a vertex-centred volume fits."""
    with uproot.open(root_path) as root_file:
        histogram = root_file["XYPseg"]
        x_edges = np.asarray(histogram.axes[0].edges(), dtype=np.float64)
        y_edges = np.asarray(histogram.axes[1].edges(), dtype=np.float64)
    try:
        crop_bounds(x_edges, truth.xpos_um, size)
        crop_bounds(y_edges, truth.ypos_um, size)
    except ValueError as error:
        raise CropSupportError(
            f"event {truth.event_id} ({root_path}): {error}"
        ) from error


def extract_volume(
    values_xyz: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    *,
    center_x_um: float,
    center_y_um: float,
    size: int = DEFAULT_VOLUME_SIZE,
    plates: int = DEFAULT_PLATES,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[int, int]]:
    """Extract [z,y,x], reproducing crop_bw.C z-bin and xy-bin conventions."""
    if values_xyz.ndim != 3:
        raise ValueError(f"expected a three-dimensional TH3 array, got {values_xyz.shape}")
    if values_xyz.shape[2] < plates + 1:
        raise ValueError(f"need ROOT z bins 2..{plates + 1}, got only {values_xyz.shape[2]}")
    x_start, x_stop, x_center_bin = crop_bounds(x_edges, center_x_um, size)
    y_start, y_stop, y_center_bin = crop_bounds(y_edges, center_y_um, size)
    # crop_bw.C calls projectHist(H3, plate=1..57), which selects ROOT bins
    # plate+1=2..58. In zero-based NumPy indexing this is [1:58].
    xyz = values_xyz[x_start:x_stop, y_start:y_stop, 1 : plates + 1]
    volume = np.rint(xyz).astype(np.int32).transpose(2, 1, 0)
    if volume.shape != (plates, size, size):
        raise RuntimeError(f"unexpected extracted shape {volume.shape}")
    return (
        volume,
        np.asarray(x_edges[x_start : x_stop + 1], dtype=np.float64),
        np.asarray(y_edges[y_start : y_stop + 1], dtype=np.float64),
        (x_center_bin, y_center_bin),
    )


def extract_available_volume(
    values_xyz: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    *,
    center_x_um: float,
    center_y_um: float,
    size: int = DEFAULT_VOLUME_SIZE,
    plates: int = DEFAULT_PLATES,
    minimum_xy_bins: int = 20,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the available part of a nominal fixed-size scan crop.

    The requested crop uses the same ROOT bin convention as ``extract_volume``:
    it is ``size`` bins around the nominal cell centre.  Where that crop exceeds
    a clipped ROOT histogram, only the physical intersection is retained.
    """
    if values_xyz.ndim != 3:
        raise ValueError(f"expected a three-dimensional TH3 array, got {values_xyz.shape}")
    if values_xyz.shape[2] < plates + 1:
        raise ValueError(f"need ROOT z bins 2..{plates + 1}, got only {values_xyz.shape[2]}")
    if values_xyz.shape[:2] != (len(x_edges) - 1, len(y_edges) - 1):
        raise ValueError("XYPseg values and XY axis edges disagree")
    if size <= 0 or size % 2:
        raise ValueError("crop size must be a positive even integer")

    def clipped_bounds(edges: np.ndarray, center: float) -> tuple[int, int]:
        bin_width = float(np.median(np.diff(edges)))
        half_span = size * bin_width / 2.0
        if center + half_span <= edges[0] or center - half_span >= edges[-1]:
            raise CropSupportError(
                f"crop of size {size} around {center} has no overlap with histogram support"
            )
        root_center_bin = int(np.searchsorted(edges, center, side="right"))
        requested_start = root_center_bin - size // 2 - 1
        requested_stop = requested_start + size
        return max(0, requested_start), min(len(edges) - 1, requested_stop)

    x_start, x_stop = clipped_bounds(x_edges, center_x_um)
    y_start, y_stop = clipped_bounds(y_edges, center_y_um)
    width = x_stop - x_start
    height = y_stop - y_start
    if width < minimum_xy_bins or height < minimum_xy_bins:
        raise CropSupportError(
            f"available part of crop is only {width}x{height} bins; need at least "
            f"{minimum_xy_bins}x{minimum_xy_bins} for a CNN crop"
        )
    volume = np.rint(
        values_xyz[x_start:x_stop, y_start:y_stop, 1 : plates + 1]
    ).astype(np.int32).transpose(2, 1, 0)
    return (
        volume,
        np.asarray(x_edges[x_start : x_stop + 1], dtype=np.float64),
        np.asarray(y_edges[y_start : y_stop + 1], dtype=np.float64),
    )


def coordinate_to_local_bin(edges: np.ndarray, value: float) -> float:
    centers = 0.5 * (edges[:-1] + edges[1:])
    return float(np.interp(value, centers, np.arange(centers.size, dtype=np.float64)))


def poisson_background_per_layer(xyseg: np.ndarray, plates: int = DEFAULT_PLATES) -> float:
    """Reproduce crop_bw.C raw-spectrum MPV and its per-layer normalization."""
    positive = np.asarray(xyseg, dtype=np.float64)
    positive = positive[positive > 50]
    if positive.size == 0:
        raise ValueError("XYseg has no bins with content > 50")
    counts, edges = np.histogram(positive, bins=200, range=(0.0, 1000.0))
    raw_mpv = float(edges[int(np.argmax(counts))])
    background_mu = raw_mpv / plates
    if background_mu <= 0:
        raise ValueError(f"invalid background mean {background_mu}")
    return background_mu


def verify_reference_crop(
    reference_path: Path,
    event_id: int,
    values_xyz: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    truth: SignalTruth,
) -> dict[str, int | bool]:
    center_x, center_y = truth.propagated_center_um()
    reconstructed, _, _, _ = extract_volume(
        values_xyz, x_edges, y_edges,
        center_x_um=center_x, center_y_um=center_y, size=32, plates=DEFAULT_PLATES,
    )
    with uproot.open(reference_path) as root_file:
        tree = root_file["samples"]
        event_ids = tree["signal_event_id"].array(library="np")
        matches = np.flatnonzero(event_ids == event_id)
        if matches.size != 1:
            raise ValueError(f"reference contains {matches.size} rows for event {event_id}")
        entry = int(matches[0])
        stored = tree["counts"].array(entry_start=entry, entry_stop=entry + 1, library="np")[0]
    difference = np.abs(reconstructed.astype(np.int64) - stored.astype(np.int64))
    return {
        "equal": bool(np.array_equal(reconstructed, stored)),
        "different_voxels": int(np.count_nonzero(difference)),
        "absolute_difference_sum": int(difference.sum()),
        "maximum_absolute_difference": int(difference.max(initial=0)),
    }


def write_hdf5(
    input_roots: Sequence[str],
    truth_by_event: dict[int, SignalTruth],
    output: Path,
    *,
    reference_samples: Path | None,
) -> dict[str, object]:
    output.parent.mkdir(parents=True, exist_ok=True)
    n = len(input_roots)
    string_dtype = h5py.string_dtype(encoding="utf-8")
    report: dict[str, object] = {"events": [], "output": str(output)}
    with h5py.File(output, "w") as hdf5:
        volumes = hdf5.create_dataset(
            "volumes_raw", shape=(n, DEFAULT_PLATES, DEFAULT_VOLUME_SIZE, DEFAULT_VOLUME_SIZE),
            dtype=np.int32, chunks=(1, 1, DEFAULT_VOLUME_SIZE, DEFAULT_VOLUME_SIZE),
            compression="gzip", compression_opts=4, shuffle=True,
        )
        event_ids = hdf5.create_dataset("event_id", shape=(n,), dtype=np.int64)
        background_ds = hdf5.create_dataset("background_mu", shape=(n,), dtype=np.float32)
        source_files = hdf5.create_dataset("source_file", shape=(n,), dtype=string_dtype)
        x_edges_ds = hdf5.create_dataset("x_edges_um", shape=(n, DEFAULT_VOLUME_SIZE + 1), dtype=np.float64)
        y_edges_ds = hdf5.create_dataset("y_edges_um", shape=(n, DEFAULT_VOLUME_SIZE + 1), dtype=np.float64)
        truth_group = hdf5.create_group("truth")
        truth_fields = {
            name: truth_group.create_dataset(name, shape=(n,), dtype=dtype)
            for name, dtype in {
                "xpos_um": np.float64, "ypos_um": np.float64,
                "projx_um": np.float64, "projy_um": np.float64,
                "ltx_mrad_raw": np.float32, "lty_mrad_raw": np.float32,
                "ltx_mrad_corrected": np.float32, "lty_mrad_corrected": np.float32,
                "energy_gev": np.float32, "plate": np.int16,
                "propagated_center_x_um": np.float64, "propagated_center_y_um": np.float64,
                "propagated_center_x_bin": np.float32, "propagated_center_y_bin": np.float32,
            }.items()
        }
        slope_bins_ds = truth_group.create_dataset("slope_xy_bins_per_z", shape=(n, 2), dtype=np.float32)

        hdf5.attrs["volume_axis_order"] = "event,z,y,x"
        hdf5.attrs["source_histogram"] = "XYPseg"
        hdf5.attrs["source_histogram_axis_order"] = "x,y,plate"
        hdf5.attrs["source_root_z_bins"] = "2..58 inclusive"
        hdf5.attrs["slope_correction"] = SLOPE_CORRECTION
        hdf5.attrs["slope_correction_formula"] = "lt_corrected = lt_raw * 1315 / 1350"
        hdf5.attrs["crop_plate_formula"] = "crop_plate = truth_plate - 3"
        hdf5.attrs["xy_bin_size_nominal_um"] = DEFAULT_XY_BIN_SIZE_UM
        hdf5.attrs["z_bin_size_um"] = DEFAULT_Z_BIN_SIZE_UM

        for index, root_path in enumerate(input_roots):
            event_id = event_id_from_signal_filename(root_path)
            if event_id not in truth_by_event:
                raise ValueError(f"event {event_id} from {root_path} absent from truth table")
            truth = truth_by_event[event_id]
            with uproot.open(root_path) as root_file:
                histogram = root_file["XYPseg"]
                values = histogram.values(flow=False)
                x_edges = np.asarray(histogram.axes[0].edges(), dtype=np.float64)
                y_edges = np.asarray(histogram.axes[1].edges(), dtype=np.float64)
                xyseg = root_file["XYseg"].values(flow=False)
            if not np.allclose(values, np.rint(values), rtol=0, atol=1e-5):
                raise ValueError(f"non-integer XYPseg contents in {root_path}")
            volume, volume_x_edges, volume_y_edges, center_bins = extract_volume(
                values, x_edges, y_edges,
                center_x_um=truth.xpos_um, center_y_um=truth.ypos_um,
            )
            # Near detector boundaries XYPseg can be clipped (and therefore its
            # axis midpoint is not the vertex).  extract_volume already checks
            # the condition that matters: a full 200x200 crop fits around the
            # truth vertex.
            axis_center = (0.5 * (x_edges[0] + x_edges[-1]), 0.5 * (y_edges[0] + y_edges[-1]))
            propagated_x, propagated_y = truth.propagated_center_um()
            reference = None
            if reference_samples is not None:
                reference = verify_reference_crop(
                    reference_samples, event_id, values, x_edges, y_edges, truth
                )
                if not reference["equal"]:
                    raise ValueError(f"event {event_id}: reference crop mismatch: {reference}")

            volumes[index] = volume
            event_ids[index] = event_id
            background_mu = poisson_background_per_layer(xyseg)
            background_ds[index] = background_mu
            source_files[index] = root_path
            x_edges_ds[index] = volume_x_edges
            y_edges_ds[index] = volume_y_edges
            values_by_name = {
                "xpos_um": truth.xpos_um, "ypos_um": truth.ypos_um,
                "projx_um": truth.projx_um, "projy_um": truth.projy_um,
                "ltx_mrad_raw": truth.ltx_mrad_raw, "lty_mrad_raw": truth.lty_mrad_raw,
                "ltx_mrad_corrected": truth.ltx_mrad_corrected,
                "lty_mrad_corrected": truth.lty_mrad_corrected,
                "energy_gev": truth.energy_gev, "plate": truth.plate,
                "propagated_center_x_um": propagated_x,
                "propagated_center_y_um": propagated_y,
                "propagated_center_x_bin": coordinate_to_local_bin(volume_x_edges, propagated_x),
                "propagated_center_y_bin": coordinate_to_local_bin(volume_y_edges, propagated_y),
            }
            for name, value in values_by_name.items():
                truth_fields[name][index] = value
            slope_bins_ds[index] = (
                truth.ltx_mrad_corrected * DEFAULT_Z_BIN_SIZE_UM / (1000.0 * DEFAULT_XY_BIN_SIZE_UM),
                truth.lty_mrad_corrected * DEFAULT_Z_BIN_SIZE_UM / (1000.0 * DEFAULT_XY_BIN_SIZE_UM),
            )
            event_report = {
                "event_id": event_id,
                "source_file": root_path,
                "source_shape_xyz": list(values.shape),
                "volume_shape_zyx": list(volume.shape),
                "volume_sum": int(volume.sum()),
                "background_mu_per_layer": background_mu,
                "vertex_root_bins": list(center_bins),
                "source_axis_center_xy_um": list(axis_center),
                "vertex_offset_from_axis_center_xy_um": [
                    truth.xpos_um - axis_center[0], truth.ypos_um - axis_center[1],
                ],
                "propagated_center_xy_um": [propagated_x, propagated_y],
                "propagated_center_xy_bin": [
                    values_by_name["propagated_center_x_bin"],
                    values_by_name["propagated_center_y_bin"],
                ],
                "reference_crop": reference,
            }
            report["events"].append(event_report)  # type: ignore[union-attr]
    report["size_bytes"] = output.stat().st_size
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", action="append", required=True, help="Repeat for each signal ROOT")
    parser.add_argument("--truth-txt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-samples-root", type=Path)
    args = parser.parse_args()
    report = write_hdf5(
        args.input_root, read_truth_txt(args.truth_txt), args.output,
        reference_samples=args.verify_samples_root,
    )
    report_path = args.output.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
