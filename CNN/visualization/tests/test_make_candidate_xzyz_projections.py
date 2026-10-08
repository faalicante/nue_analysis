import numpy as np

from visualization.make_candidate_xzyz_projections import visual_qhot_projections


def test_visual_qhot_projections_uses_requested_poisson_tail_probability() -> None:
    raw = np.zeros((57, 8, 8), dtype=float)
    raw[10, 4, 4] = 7.0

    projections, threshold = visual_qhot_projections(
        raw, 1.0, 2, 2, 4, smooth=False, poisson_tail_probability=0.1
    )

    assert threshold == 3
    assert projections[0].sum() == 6.0
    assert projections[1].sum() == 6.0


def test_visual_qhot_projections_smooths_full_volume_before_crop() -> None:
    raw = np.zeros((57, 8, 8), dtype=float)
    raw[10, 2, 1] = 100.0

    raw_projections, _ = visual_qhot_projections(
        raw, 0.1, 2, 2, 4, smooth=False, poisson_tail_probability=0.5
    )
    smooth_projections, _ = visual_qhot_projections(
        raw, 0.1, 2, 2, 4, smooth=True, poisson_tail_probability=0.5
    )

    assert raw_projections[0].sum() == 0.0
    assert smooth_projections[0].sum() > 0.0
