from __future__ import annotations

import numpy as np
import pytest
import torch

from training.cnn21d.losses import MultiTaskLoss
from training.cnn21d.metrics import slopes_to_angles
from training.cnn21d.model import TwoHeadCNN21D
from training.train_cnn21d import build_task_targets, sample_types_to_labels


def test_shapes_through_every_layer() -> None:
    model = TwoHeadCNN21D()
    features, shapes = model.extract_features(torch.randn(2, 1, 57, 20, 20), return_shapes=True)
    assert shapes == {
        "input": (2, 1, 57, 20, 20),
        "block1": (2, 16, 57, 10, 10),
        "block2": (2, 32, 57, 5, 5),
        "block3": (2, 64, 28, 5, 5),
        "block4": (2, 128, 14, 5, 5),
        "global_features": (2, 256),
    }
    assert features.shape == (2, 256)
    outputs = model(torch.randn(2, 1, 57, 20, 20))
    assert outputs["presence_logit"].shape == (2,)
    assert outputs["slope_xy"].shape == (2, 2)


def test_forward_backward_batch() -> None:
    model = TwoHeadCNN21D()
    criterion = MultiTaskLoss(presence_weight=0.7, slope_weight=1.3)
    outputs = model(torch.randn(4, 1, 57, 20, 20))
    loss = criterion(outputs, torch.tensor([0.0, 1.0, 0.0, 1.0]), torch.randn(4, 2))["total"]
    loss.backward()
    assert torch.isfinite(loss)
    assert all(parameter.grad is not None for parameter in model.parameters())


def test_regression_ignores_negative_targets_and_gradients() -> None:
    criterion = MultiTaskLoss()
    prediction = torch.randn(4, 2, requires_grad=True)
    outputs = {"presence_logit": torch.zeros(4, requires_grad=True), "slope_xy": prediction}
    presence = torch.tensor([1.0, 0.0, 1.0, 0.0])
    targets = torch.randn(4, 2)
    first = criterion(outputs, presence, targets)["slope"]
    changed = targets.clone(); changed[1] = 1.0e6; changed[3] = -1.0e6
    second = criterion(outputs, presence, changed)["slope"]
    torch.testing.assert_close(first, second, rtol=0, atol=0)
    first.backward()
    assert torch.count_nonzero(prediction.grad[[1, 3]]) == 0
    assert torch.count_nonzero(prediction.grad[[0, 2]]) > 0


def test_batch_without_positives_has_zero_safe_regression_loss() -> None:
    criterion = MultiTaskLoss()
    prediction = torch.randn(3, 2, requires_grad=True)
    outputs = {"presence_logit": torch.randn(3, requires_grad=True), "slope_xy": prediction}
    losses = criterion(outputs, torch.zeros(3), torch.full((3, 2), 999.0))
    assert losses["slope"].item() == 0.0
    losses["total"].backward()
    torch.testing.assert_close(prediction.grad, torch.zeros_like(prediction))


def test_theta_balanced_regression_uses_configured_weights() -> None:
    criterion = MultiTaskLoss(
        regression_theta_bin_edges_mrad=[0.0, 10.0, 100.0],
        regression_theta_bin_weights=[1.0, 3.0],
    )
    prediction = torch.zeros(2, 2, requires_grad=True)
    target = torch.tensor([[0.1, 0.0], [1.0, 0.0]])
    losses = criterion(
        {"presence_logit": torch.zeros(2), "slope_xy": prediction},
        torch.ones(2),
        target,
    )
    # Per-example component means are 0.0025 and 0.25. Weight-normalized
    # Smooth-L1 is therefore (1*0.0025 + 3*0.25)/(1+3).
    torch.testing.assert_close(losses["slope"], torch.tensor(0.188125))
    torch.testing.assert_close(
        losses["effective_regression_weight"], torch.tensor(4.0)
    )


def test_low_angle_signal_is_ignored_by_bce_but_used_by_regression() -> None:
    sample_type = torch.tensor([1, 2, 2, 2])
    # theta is approximately 0, 10, exactly 15, and 20 mrad respectively.
    slope = torch.tensor([
        [0.0, 0.0], [0.270009, 0.0],
        [27.0 * torch.tan(torch.tensor(0.015)).item(), 0.0], [0.540072, 0.0],
    ])
    labels, presence_mask, regression_mask = build_task_targets(
        sample_type, slope, positive_sample_types=(2,),
        presence_positive_theta_min_mrad=15.0,
    )
    torch.testing.assert_close(labels, torch.tensor([0.0, 0.0, 0.0, 1.0]))
    assert presence_mask.tolist() == [True, False, False, True]
    assert regression_mask.tolist() == [False, True, True, True]

    logits = torch.zeros(4, requires_grad=True)
    prediction = torch.zeros(4, 2, requires_grad=True)
    criterion = MultiTaskLoss()
    losses = criterion(
        {"presence_logit": logits, "slope_xy": prediction}, labels, slope,
        presence_mask=presence_mask, regression_mask=regression_mask,
    )
    losses["total"].backward()
    assert logits.grad[1].item() == 0.0
    assert logits.grad[2].item() == 0.0
    assert torch.count_nonzero(prediction.grad[0]).item() == 0
    assert torch.count_nonzero(prediction.grad[1]).item() > 0
    assert torch.count_nonzero(prediction.grad[2]).item() > 0
    assert torch.count_nonzero(prediction.grad[3]).item() > 0


def test_low_angle_signal_can_be_bce_negative_and_regression_positive() -> None:
    sample_type = torch.tensor([1, 2, 2, 2])
    slope = torch.tensor([
        [0.0, 0.0], [0.05, 0.0], [0.27, 0.0], [0.54, 0.0],
    ])
    labels, presence_mask, regression_mask = build_task_targets(
        sample_type,
        slope,
        positive_sample_types=(2,),
        presence_positive_theta_min_mrad=10.0,
        presence_low_angle_policy="negative",
    )
    # Hard and the two signal examples at/below 10 mrad are BCE negatives;
    # every signal example, including low-angle signal, remains regression-valid.
    torch.testing.assert_close(labels, torch.tensor([0.0, 0.0, 0.0, 1.0]))
    assert presence_mask.tolist() == [True, True, True, True]
    assert regression_mask.tolist() == [False, True, True, True]

    logits = torch.zeros(4, requires_grad=True)
    prediction = torch.zeros(4, 2, requires_grad=True)
    losses = MultiTaskLoss()(
        {"presence_logit": logits, "slope_xy": prediction},
        labels,
        slope,
        presence_mask=presence_mask,
        regression_mask=regression_mask,
    )
    losses["total"].backward()
    assert torch.count_nonzero(logits.grad).item() == 4
    assert torch.count_nonzero(prediction.grad[0]).item() == 0
    assert torch.count_nonzero(prediction.grad[1:]).item() > 0


def test_transition_band_is_ignored_by_bce_but_kept_for_regression() -> None:
    sample_type = torch.tensor([1, 2, 2, 2, 2])
    # Signal examples lie below 5 mrad, in 5-10 mrad, at 10 mrad, and above 10 mrad.
    slope = torch.tensor([
        [0.0, 0.0], [0.10, 0.0], [0.20, 0.0], [0.27, 0.0], [0.40, 0.0],
    ])
    labels, presence_mask, regression_mask = build_task_targets(
        sample_type,
        slope,
        positive_sample_types=(2,),
        presence_positive_theta_min_mrad=10.0,
        presence_low_angle_policy="negative",
        presence_negative_theta_max_mrad=5.0,
    )
    torch.testing.assert_close(labels, torch.tensor([0.0, 0.0, 0.0, 0.0, 1.0]))
    assert presence_mask.tolist() == [True, True, False, False, True]
    assert regression_mask.tolist() == [False, True, True, True, True]


def test_invalid_low_angle_policy_is_rejected() -> None:
    with pytest.raises(ValueError, match="presence_low_angle_policy"):
        build_task_targets(
            torch.tensor([2]), torch.tensor([[0.1, 0.0]]), (2,), 10.0, "bad"
        )


def test_angle_and_axis_convention() -> None:
    theta, phi = slopes_to_angles(np.asarray([[27.0, 0.0], [0.0, 27.0]]))
    np.testing.assert_allclose(theta, 1000.0 * np.pi / 4.0)
    np.testing.assert_allclose(phi, [0.0, np.pi / 2.0])


def test_binary_labels_come_from_sample_type_and_exclude_poisson_by_config() -> None:
    sample_type = torch.tensor([0, 1, 2, 1, 2])
    labels = sample_types_to_labels(sample_type, positive_sample_types=(2,))
    torch.testing.assert_close(labels, torch.tensor([0.0, 0.0, 1.0, 0.0, 1.0]))
    include_sample_types = (1, 2)
    eligible = torch.zeros_like(sample_type, dtype=torch.bool)
    for value in include_sample_types:
        eligible |= sample_type == value
    assert eligible.tolist() == [False, True, True, True, True]
