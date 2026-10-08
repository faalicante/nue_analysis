from __future__ import annotations

import numpy as np

from studies.evaluate_peripheral_background import distance_from_window_to_trajectory
from scanning.scan_cnn21d_volumes import (
    augment_volume_background,
    cluster_contains_truth,
    connected_components,
    sliding_positions,
    truth_reference_bins,
)


def test_scan_poisson_background_preserves_zero_voxels() -> None:
    volume = np.full((3, 4, 5), 4, dtype=np.int16)
    volume[:, 1, 2] = 0
    augmented, effective_mu, added_mu = augment_volume_background(
        volume, 4.0, 7.5, seed=17
    )
    assert effective_mu == 7.5
    assert added_mu == 3.5
    np.testing.assert_array_equal(augmented[volume == 0], 0.0)
    assert np.any(augmented[volume > 0] > volume[volume > 0])
    repeated, _, _ = augment_volume_background(volume, 4.0, 7.5, seed=17)
    np.testing.assert_array_equal(augmented, repeated)


def test_window_distance_from_trajectory() -> None:
    trajectory_x = np.asarray([5.0, 15.0, 25.0])
    trajectory_y = np.asarray([5.0, 5.0, 5.0])
    assert distance_from_window_to_trajectory(10, 0, 10, trajectory_x, trajectory_y) == 0.0
    assert np.isclose(
        distance_from_window_to_trajectory(10, 20, 10, trajectory_x, trajectory_y),
        15.0,
    )


def test_stride_grid_for_200_by_200_volume() -> None:
    positions = sliding_positions(200, 20, 10)
    assert positions == list(range(0, 181, 10))
    assert len(positions) == 19


def test_stride_grid_includes_last_edge_when_not_divisible() -> None:
    assert sliding_positions(25, 10, 6) == [0, 6, 12, 15]


def test_eight_connected_clustering_and_truth_matching() -> None:
    mask = np.zeros((4, 4), dtype=bool)
    mask[0, 0] = True
    mask[1, 1] = True
    mask[3, 3] = True
    components = connected_components(mask)
    assert sorted(map(len, components)) == [1, 2]
    diagonal_component = next(component for component in components if len(component) == 2)
    positions = [0, 10, 20, 30]
    assert cluster_contains_truth(diagonal_component, 15.0, 15.0, positions, positions, 20)
    assert not cluster_contains_truth(diagonal_component, 39.0, 39.0, positions, positions, 20)


def test_truth_reference_bins_can_use_24_plate_limit() -> None:
    class FakeDataset:
        def __init__(self, values: list[float]) -> None:
            self.values = values

        def __getitem__(self, index: int) -> float:
            return self.values[index]

    class FakeSource:
        def __getitem__(self, key: str) -> object:
            assert key == "truth"
            return {
                "plate": FakeDataset([11]),
                "xpos_um": FakeDataset([100.0]),
                "ypos_um": FakeDataset([200.0]),
                "ltx_mrad_corrected": FakeDataset([10.0]),
                "lty_mrad_corrected": FakeDataset([-5.0]),
            }

    edges = np.arange(-5.0, 406.0, 10.0)
    x_bin, y_bin = truth_reference_bins(
        FakeSource(), 0, edges, edges, plates=57, propagation_limit=24  # type: ignore[arg-type]
    )
    assert np.isclose(x_bin, 26.2)
    assert np.isclose(y_bin, 11.9)
