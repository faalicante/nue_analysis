# Model and training

Run Python tools from the `CNN` project root as `python -m training.<script_name>` (omit `.py`). Run shell launchers with `bash training/<script_name>.sh`.

| Code | Purpose |
|---|---|
| [cnn_dataset.py](cnn_dataset.py) | Load HDF5 samples, normalize inputs, apply coherent crop/slope augmentation and build DataLoaders. |
| [train_cnn21d.py](train_cnn21d.py) | Train and evaluate the two-headed CNN, saving classification and regression checkpoints. |
| [train_gap5_10_mu_regimes.sh](train_gap5_10_mu_regimes.sh) | Launch low- and high-background-mu training experiments. |
| [train_gap5_10_mu_aug_global.sh](train_gap5_10_mu_aug_global.sh) | Launch training with global Poisson-background augmentation. |
| [train_fit_theta_normalization_tests.sh](train_fit_theta_normalization_tests.sh) | Compare three normalization modes on the fit-selected hard dataset. |

`cnn21d/model.py` defines the network, `cnn21d/losses.py` the multi-task loss, and `cnn21d/metrics.py` the evaluation metrics. `cnn21d/__init__.py` exposes the shared model/loss API. `configs/` contains the experiment settings. See [training](CNN21D_TRAINING.md) and [input pipeline](CNN_DATALOADER.md).
