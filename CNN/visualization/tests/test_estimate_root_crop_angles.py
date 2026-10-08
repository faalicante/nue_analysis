from unittest.mock import patch

import numpy as np

from studies.estimate_root_crop_angles import fit_record, sliding_positions, theta_mrad


def test_cnn_windows_cover_a_32_bin_crop():
    assert sliding_positions(32, 20, 4) == [0, 4, 8, 12]


def test_theta_conversion_uses_the_xypseg_aspect_ratio():
    assert np.isclose(theta_mrad(np.asarray([0.0, 0.0])), 0.0)
    assert np.isclose(theta_mrad(np.asarray([0.0, 0.27])), 10.0, atol=0.01)


def test_fit_record_keeps_crop_when_no_qhot_component_is_found():
    missing_fit = {
        "start_layer": None,
        "start_x_um": None,
        "start_y_um": None,
        "start_component_voxels": 0,
        "start_component_layers": 0,
        "start_quality": "not_found",
        "poisson_hot_count_threshold": 13,
        "component_selection": "not_found",
    }
    with patch("studies.estimate_root_crop_angles.estimate_shower_start", return_value=missing_fit):
        row = fit_record(
            np.zeros((57, 32, 32)), 4.0, np.asarray([0.1, 0.2]), (16.0, 16.0)
        )

    assert row["fit_quality"] == "not_found"
    assert row["fit_theta_mrad"] == ""
    assert row["fit_slope_x_bins_per_layer"] == ""
