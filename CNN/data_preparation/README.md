# Data preparation

Run Python tools from the `CNN` project root as `python -m data_preparation.<script_name>` (omit `.py`). Run shell launchers with `bash data_preparation/<script_name>.sh`.

| Code | Purpose |
|---|---|
| [crop_bw.C](crop_bw.C) | Export raw ROOT crops for signal, hard background and Poisson background. |
| [loop.py](loop.py) | Legacy MC event launcher and zScale collector; currently limited to three selected samples. |
| [root_to_hdf5.py](root_to_hdf5.py) | Convert raw ROOT source-region TTrees to a validated HDF5 dataset. |
| [extract_signal_crops_from_scan_shards.py](extract_signal_crops_from_scan_shards.py) | Rebuild centred signal ROOT crops directly from 57x200x200 scan shards. |
| [align_hdf5_splits_to_reference.py](align_hdf5_splits_to_reference.py) | Align HDF5 train/validation/test groups to an existing manifest. |
| [build_scan_split_hdf5.py](build_scan_split_hdf5.py) | Collect sharded scan volumes into original event-group validation/test splits. |
| [build_fit_theta_hard_dataset.py](build_fit_theta_hard_dataset.py) | Build a training HDF5 with hard crops selected by coherent Poisson-hot fit. |
| [cache_local_mpv_z_stats.py](cache_local_mpv_z_stats.py) | Cache per-layer local MPV and lower-side width for a CNN training HDF5. |

Build the ROOT producer from the project root with `make`, then run `./data_preparation/run.exe`. ROOT and its development tools are required. See [conversion details](ROOT_TO_HDF5.md) and [producer details](SOURCE_REGION_EXPORT.md).
