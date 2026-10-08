from __future__ import annotations

import pytest

from scanning.produce_signal_scan_shards import (
    chunked,
    loop_selection_failures,
    render_root_path,
    theta_gt5_selection_failures,
    training_selection_failures,
)
from scanning.root_xypseg_to_scan_hdf5 import SignalTruth


def test_render_signal_root_path() -> None:
    template = "/eos/data/b000021.0.0.{event_plus_one}.trk.root"
    assert render_root_path(template, 25).endswith(".26.trk.root")
    assert render_root_path("event_{event}.root", 25) == "event_25.root"


def test_chunked_events() -> None:
    assert chunked([1, 2, 3, 4, 5], 2) == [[1, 2], [3, 4], [5]]
    with pytest.raises(ValueError):
        chunked([1], 0)


def make_truth(*, energy: float = 100.0, plate: int = 20, tx: float = 10.0) -> SignalTruth:
    return SignalTruth(1, 0, 0, 0, 0, tx, 0, energy, plate)


def test_loop_selection_reproduces_strict_cuts() -> None:
    assert loop_selection_failures(make_truth()) == []
    assert "energy_not_gt_80_gev" in loop_selection_failures(make_truth(energy=80.0))
    assert "crop_plate_not_lt_40" in loop_selection_failures(make_truth(plate=43))
    assert "theta_not_gt_5_mrad" in loop_selection_failures(make_truth(tx=5.0))
    assert "theta_not_lt_25_mrad" in loop_selection_failures(make_truth(tx=26.0))


def test_default_theta_selection_ignores_energy_plate_and_has_no_upper_cut() -> None:
    assert theta_gt5_selection_failures(make_truth(energy=1.0, plate=60, tx=100.0)) == []
    assert theta_gt5_selection_failures(make_truth(tx=5.0)) == ["theta_not_gt_5_mrad"]


def test_training_selection_has_energy_plate_and_lower_theta_but_no_upper_theta() -> None:
    assert training_selection_failures(make_truth(tx=100.0)) == []
    assert "energy_not_gt_80_gev" in training_selection_failures(make_truth(energy=80.0))
    assert "crop_plate_not_lt_40" in training_selection_failures(make_truth(plate=43))
    assert "theta_not_gt_5_mrad" in training_selection_failures(make_truth(tx=5.0))
