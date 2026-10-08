from __future__ import annotations

import numpy as np

from scanning.root_xypseg_to_scan_hdf5 import (
    SLOPE_CORRECTION,
    SignalTruth,
    crop_bounds,
    event_id_from_signal_filename,
    extract_available_volume,
    extract_volume,
    poisson_background_per_layer,
)


def test_signal_filename_event_mapping() -> None:
    assert event_id_from_signal_filename("b000021.0.0.26.trk.root") == 25
    assert event_id_from_signal_filename("root://host/path/b000021.0.0.6009.trk.root") == 6008


def test_slope_correction_and_propagated_center() -> None:
    truth = SignalTruth(25, 242024.64, 72013.84, 0, 0, 13.82, -0.77, 935.33, 15)
    assert np.isclose(truth.ltx_mrad_corrected, 13.82 * SLOPE_CORRECTION)
    assert np.allclose(
        truth.propagated_center_um(),
        (242024.64 + 13.82 * SLOPE_CORRECTION * 27.0,
         72013.84 - 0.77 * SLOPE_CORRECTION * 27.0),
    )

    late_truth = SignalTruth(6008, 364657.26, 50673.81, 0, 0, -16.76, -44.48, 88.38, 26)
    factor = 0.5 * 1350.0 * 34 / 1000.0  # crop plate = 26 - 3 = 23
    assert np.allclose(
        late_truth.propagated_center_um(),
        (364657.26 - 16.76 * SLOPE_CORRECTION * factor,
         50673.81 - 44.48 * SLOPE_CORRECTION * factor),
    )


def test_even_crop_matches_root_center_convention_and_z_offset() -> None:
    values = np.arange(10 * 10 * 60, dtype=np.float32).reshape(10, 10, 60)
    edges = np.arange(11, dtype=np.float64)
    start, stop, center_bin = crop_bounds(edges, 5.2, 4)
    assert (start, stop, center_bin) == (3, 7, 6)
    volume, x_edges, y_edges, bins = extract_volume(
        values, edges, edges, center_x_um=5.2, center_y_um=5.2, size=4, plates=3
    )
    expected = values[3:7, 3:7, 1:4].transpose(2, 1, 0).astype(np.int32)
    np.testing.assert_array_equal(volume, expected)
    np.testing.assert_array_equal(x_edges, edges[3:8])
    np.testing.assert_array_equal(y_edges, edges[3:8])
    assert bins == (6, 6)


def test_extract_available_volume_preserves_a_clipped_xy_area() -> None:
    values = np.arange(24 * 96 * 60, dtype=np.float32).reshape(24, 96, 60)
    x_edges = np.arange(25, dtype=np.float64) * 50.0
    y_edges = np.arange(97, dtype=np.float64) * 50.0
    volume, local_x_edges, local_y_edges = extract_available_volume(
        values,
        x_edges,
        y_edges,
        center_x_um=600.0,
        center_y_um=2_400.0,
        plates=3,
        minimum_xy_bins=20,
    )
    assert volume.shape == (3, 96, 24)
    np.testing.assert_array_equal(volume, values[:, :, 1:4].transpose(2, 1, 0))
    np.testing.assert_array_equal(local_x_edges, x_edges)
    np.testing.assert_array_equal(local_y_edges, y_edges)


def test_extract_available_volume_is_limited_to_nominal_200_bin_crop() -> None:
    values = np.arange(240 * 230 * 5, dtype=np.float32).reshape(240, 230, 5)
    x_edges = np.arange(241, dtype=np.float64) * 50.0
    y_edges = np.arange(231, dtype=np.float64) * 50.0

    volume, local_x_edges, local_y_edges = extract_available_volume(
        values,
        x_edges,
        y_edges,
        center_x_um=6_000.0,
        center_y_um=5_750.0,
        plates=3,
        minimum_xy_bins=20,
    )

    assert volume.shape == (3, 200, 200)
    np.testing.assert_array_equal(volume, values[20:220, 15:215, 1:4].transpose(2, 1, 0))
    np.testing.assert_array_equal(local_x_edges, x_edges[20:221])
    np.testing.assert_array_equal(local_y_edges, y_edges[15:216])


def test_poisson_background_matches_root_histogram_mpv_convention() -> None:
    xyseg = np.concatenate((np.full(20, 197.0), np.full(5, 212.0), [0.0, 1500.0]))
    # 197 falls in the ROOT-style [195, 200) spectrum bin.
    assert poisson_background_per_layer(xyseg, plates=57) == 195.0 / 57.0
