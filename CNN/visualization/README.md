# Visualization, validation and utilities

Run Python tools from the `CNN` project root as `python -m visualization.<script_name>` (omit `.py`). Run shell launchers with `bash visualization/<script_name>.sh`.

| Code | Purpose |
|---|---|
| [make_background_candidate_gifs.py](make_background_candidate_gifs.py) | Create one smoothed 57-layer GIF for each background scan candidate. |
| [make_signal_scan_gifs.py](make_signal_scan_gifs.py) | Create smoothed GIFs for missed and clearly detected signal scan events. |
| [make_original_crop_gifs.py](make_original_crop_gifs.py) | Create smoothed 57-layer GIFs directly from original 32x32 ROOT crops. |
| [make_root_crop_gifs_from_angles.py](make_root_crop_gifs_from_angles.py) | Create smoothed layer GIFs for ROOT crops selected from an angle CSV. |
| [make_candidate_xzyz_projections.py](make_candidate_xzyz_projections.py) | Export XZ/YZ q-hot projections with CNN and geometric slopes. |
| [make_root_crop_xzyz_projections.py](make_root_crop_xzyz_projections.py) | Render XZ/YZ projections of ROOT crops selected from an angle CSV. |
| [plot_candidate_xzyz.py](plot_candidate_xzyz.py) | Plot Poisson-hot XZ/YZ projections for selected scan-candidate ranks. |
| [tree_to_images.py](tree_to_images.py) | Render raw and smoothed/cropped images from a ROOT TTree counts branch. |
| [plot_background_scan_candidates.py](plot_background_scan_candidates.py) | Render 57-layer contact sheets for truth-free background scan candidates. |
| [plot_missed_signal_scan_events.py](plot_missed_signal_scan_events.py) | Render truth-centred CNN windows for signal events missed by the stride scan. |
| [plot_hard_theta_distribution.py](plot_hard_theta_distribution.py) | Plot the predicted-angle distribution for the unbiased hard control sample. |
| [plot_pilot_diagnostics.py](plot_pilot_diagnostics.py) | Generate detailed diagnostic plots from the no-Poisson pilot checkpoint. |
| [candidates.sh](candidates.sh) | Generate projections and GIFs for a configured list of data bricks and models. |
| [make_background_candidate_gifs.sh](make_background_candidate_gifs.sh) | Launch candidate GIF generation for one data brick/model. |
| [make_candidate_xzyz_projections.sh](make_candidate_xzyz_projections.sh) | Launch candidate XZ/YZ projections for one data brick/model. |
| [qa_cnn_dataset.py](qa_cnn_dataset.py) | Visual QA for central versus training-augmented 2+1D CNN inputs. |
| [smoke_test_cnn_dataloader.py](smoke_test_cnn_dataloader.py) | Load one real batch from every HDF5 split and print a JSON summary. |
| [rename_brick_files.py](rename_brick_files.py) | Rename brick identifiers in filenames; preview by default and apply only with --apply. |

`tests/` contains the automated regression tests. [Presentation sources](presentations/README.md) are separate from final slide decks in `output/`. The three shell launchers generate candidate graphics; `rename_brick_files.py` is a filename maintenance utility.
