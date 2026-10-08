import pytest

from candidates.export_vertex_search_seeds import (
    candidate_source_location,
    project_xy,
    vertex_interval,
)


def test_vertex_interval_uses_mc_residual_and_physical_constraint():
    # With visible L=24 and r16=1, r84=9, vertex is [24-3-9, 24-3-1].
    assert vertex_interval(24, 1, 9) == (12, 20, 12, 20)


def test_vertex_interval_is_clipped_at_first_detector_layer():
    assert vertex_interval(8, 1, 9) == (-4, 4, 1, 4)


def test_projection_runs_backwards_in_z():
    x, y = project_xy(1000.0, 2000.0, 10.0, -20.0, 20, 10)
    assert x == 865.0
    assert y == 2270.0


def test_signal_manifest_source_is_resolved_by_event_id():
    assert candidate_source_location(
        {"event_id": 42}, {42: ("signal_scan_test.h5", 7)}
    ) == ("signal_scan_test.h5", 7)
    with pytest.raises(ValueError, match="pass --volumes"):
        candidate_source_location({"event_id": 42}, {})
