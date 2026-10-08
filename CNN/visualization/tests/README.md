# Automated validation

Run `python -m pytest -q` from the project root. These tests use synthetic temporary data and do not train production models or modify preserved datasets.

| Test file | Coverage |
|---|---|
| [test_root_to_hdf5.py](test_root_to_hdf5.py) | ROOT conversion, count preservation, metadata, grouping and validation errors. |
| [test_cnn_dataset.py](test_cnn_dataset.py) | Normalization, crop/augmentation geometry and worker-local DataLoader handles. |
| [test_cnn21d.py](test_cnn21d.py) | Model dimensions, gradients, task targets and loss masking/weighting. |
| [test_root_xypseg_to_scan_hdf5.py](test_root_xypseg_to_scan_hdf5.py) | ROOT volume extraction, bin conventions and scan-volume metadata. |
| [test_background_scan_shards.py](test_background_scan_shards.py) | Background/data cell geometry, custom offsets, partial cells and shard writing. |
| [test_produce_signal_scan_shards.py](test_produce_signal_scan_shards.py) | Signal production selection and shard utilities. |
| [test_produce_data_scan_shards.py](test_produce_data_scan_shards.py) | Brick identifiers, selected-cell lists and data provenance. |
| [test_scan_cnn21d_volumes.py](test_scan_cnn21d_volumes.py) | Sliding windows, clustering and scan geometry. |
| [test_estimate_root_crop_angles.py](test_estimate_root_crop_angles.py) | Crop retention when no coherent geometric fit is found. |
| [test_export_vertex_search_seeds.py](test_export_vertex_search_seeds.py) | Direction conversion and vertex-search seed coordinates. |
| [test_make_background_candidate_gifs.py](test_make_background_candidate_gifs.py) | Candidate display and shower-start geometry. |
| [test_make_candidate_xzyz_projections.py](test_make_candidate_xzyz_projections.py) | Poisson-hot projection behavior. |
