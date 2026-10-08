# Volume production and CNN scanning

Run Python tools from the `CNN` project root as `python -m scanning.<script_name>` (omit `.py`). Run shell launchers with `bash scanning/<script_name>.sh`.

| Code | Purpose |
|---|---|
| [root_xypseg_to_scan_hdf5.py](root_xypseg_to_scan_hdf5.py) | Extract vertex-centred 57x200x200 scan volumes from signal XYPseg TH3 files. |
| [root_background_xypseg_to_scan_hdf5.py](root_background_xypseg_to_scan_hdf5.py) | Extract physical 57x200x200 background cells from XYPseg ROOT files. |
| [produce_signal_scan_shards.py](produce_signal_scan_shards.py) | Produce restartable HDF5 shards of 57x200x200 signal scan volumes. |
| [produce_background_scan_shards.py](produce_background_scan_shards.py) | Produce restartable HDF5 shards for full-cell background CNN scans. |
| [produce_data_scan_shards.py](produce_data_scan_shards.py) | Produce restartable 57x200x200 HDF5 shards for real-data CNN scans. |
| [scan_cnn21d_volumes.py](scan_cnn21d_volumes.py) | Scan 57x200x200 HDF5 volumes with overlapping CNN crops and cluster detections. |
| [scan_loop.sh](scan_loop.sh) | Launch the configured gap5_10_mu_high scans for a list of data bricks. |
| [scan_cnn21d_loop.sh](scan_cnn21d_loop.sh) | Scan a data brick with the selected gap/high-mu model. |
| [scan_cnn21d_loop2.sh](scan_cnn21d_loop2.sh) | Scan a data brick with the fitlt5_mu_residual model. |
| [scan_high_mu_hard_background.sh](scan_high_mu_hard_background.sh) | Produce, scan and summarize a separate hard-background ROOT production. |
| [scan_fit_theta_normalization_tests.sh](scan_fit_theta_normalization_tests.sh) | Compare signal/background scans for the three normalization modes. |
| [evaluate_mu_regime_scans.sh](evaluate_mu_regime_scans.sh) | Evaluate low/high-mu models on signal and background scan volumes. |
| [test_scan.sh](test_scan.sh) | Run signal-volume inference with the resumed dual-checkpoint model; this is not a unit test. |

Cells are 10,000 µm wide. Background grid offsets refer to the lower edges of cell (1,1): b21 defaults to (200000, 4500) µm; the separate b24 production uses (5729, 196206) µm. Data cell centers are (i×10000, j×10000) µm. Signal scan volumes are centered on MC vertices.
