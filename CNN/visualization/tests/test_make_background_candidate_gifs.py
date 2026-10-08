import numpy as np

from visualization.make_background_candidate_gifs import (
    centered_window_start,
    estimate_shower_start,
    fitted_display_crop_size,
)


def test_centered_window_start_and_boundary_clamp() -> None:
    assert centered_window_start(80, 20, 40, 200) == 70
    assert centered_window_start(0, 20, 40, 200) == 0
    assert centered_window_start(180, 20, 40, 200) == 160


def test_fitted_display_crop_size_for_partial_cells() -> None:
    assert fitted_display_crop_size(40, 20, 43, 147) == 40
    assert fitted_display_crop_size(40, 20, 38, 50) == 38


def test_estimate_shower_start_from_coherent_hot_component() -> None:
    raw = np.zeros((57, 40, 40), dtype=np.float64)
    raw[5, 10, 11] = 20
    raw[6, 10, 11] = 20
    raw[7, 10, 11] = 20
    raw[7, 10, 12] = 20
    raw[7, 11, 11] = 20
    raw[7, 11, 12] = 20
    raw[8, 10, 11] = 20
    raw[8, 10, 12] = 20
    x_centers = np.arange(40, dtype=np.float64) * 50.0 + 25.0
    y_centers = np.arange(40, dtype=np.float64) * 50.0 + 25.0

    result = estimate_shower_start(raw, 1.0, x_centers, y_centers)

    assert result["start_layer"] == 6
    assert result["start_root_z_bin"] == 7
    assert result["start_x_um"] == x_centers[11]
    assert result["start_y_um"] == y_centers[10]
    assert result["start_component_layers"] == 4
    assert result["start_component_voxels"] == 8
    assert result["start_quality"] == "coherent"


def test_estimate_shower_start_joins_gapped_ridge_matching_cnn_direction() -> None:
    raw = np.zeros((57, 40, 40), dtype=np.float64)
    # A brighter but horizontal object should not hide this fragmented diagonal
    # ridge when the CNN direction points along the latter.
    for z, x, y in ((10, 8, 7), (12, 9, 8), (14, 10, 9), (16, 11, 10)):
        raw[z, y, x] = 20
    for z in range(25, 31):
        raw[z, 25, 25] = 40
    x_centers = np.arange(40, dtype=np.float64) * 50.0 + 25.0
    y_centers = np.arange(40, dtype=np.float64) * 50.0 + 25.0

    result = estimate_shower_start(
        raw,
        1.0,
        x_centers,
        y_centers,
        predicted_slope_xy=np.asarray([0.5, 0.5]),
        candidate_center_xy_um=np.asarray([1000.0, 1000.0]),
    )

    assert result["component_selection"] == "cnn_consistent_qhot_ridge"
    assert result["start_layer"] == 11
    assert result["start_component_layers"] == 4
    np.testing.assert_allclose(
        [result["selected_component_slope_x"], result["selected_component_slope_y"]],
        [0.5, 0.5],
        atol=1.0e-12,
    )
