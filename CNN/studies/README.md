# Experimental studies and comparisons

Run Python tools from the `CNN` project root as `python -m studies.<script_name>` (omit `.py`). Run shell launchers with `bash studies/<script_name>.sh`.

| Code | Purpose |
|---|---|
| [evaluate_angle15_pilots.py](evaluate_angle15_pilots.py) | Compare angle-threshold pilot operating points and regression by true theta. |
| [evaluate_angle15_seed_stability.py](evaluate_angle15_seed_stability.py) | Evaluate angle>15 classification stability on complete validation/test splits. |
| [compare_presence_theta_thresholds.py](compare_presence_theta_thresholds.py) | Compare presence checkpoints trained with different signal theta thresholds. |
| [evaluate_intermediate_operating_point.py](evaluate_intermediate_operating_point.py) | Choose a presence threshold on validation and apply it once to test. |
| [evaluate_l5_pilot.py](evaluate_l5_pilot.py) | Compare the low-angle pilot with the current full model on held-out data. |
| [evaluate_jitter_offsets.py](evaluate_jitter_offsets.py) | Evaluate the angle15 crop20 pilot at fixed scan-grid miscenterings. |
| [evaluate_peripheral_background.py](evaluate_peripheral_background.py) | Calibrate scan scores using signal-centred maps away from the MC trajectory. |
| [analyze_shower_size_features.py](analyze_shower_size_features.py) | Compare Poisson-robust shower-size features in signal and hard samples. |
| [evaluate_scan_score_size_operating_points.py](evaluate_scan_score_size_operating_points.py) | Evaluate CNN score thresholds before and after the shower-size post-filter. |
| [estimate_candidate_yield.py](estimate_candidate_yield.py) | Estimate candidate yield by folding pilot response with cumulative angle efficiencies. |
| [summarize_scan_volume_predictions.py](summarize_scan_volume_predictions.py) | Summarize event-level scan efficiency and candidate multiplicity versus threshold. |
| [summarize_signal_angle_bins.py](summarize_signal_angle_bins.py) | Report pilot performance in signal theta bins for the preselected t5 sample. |
| [fit_signal_centroids.py](fit_signal_centroids.py) | Fit straight shower axes directly to significant voxels in the raw HDF5 volumes. |
| [estimate_root_crop_angles.py](estimate_root_crop_angles.py) | Estimate CNN and Poisson-hot fit angles in original 57x32x32 ROOT crops. |
| [summarize_root_crop_angles.py](summarize_root_crop_angles.py) | Summarize and plot CNN/Poisson-fit angles from ROOT crop angle CSV files. |
| [compare_event938_crop_centers.py](compare_event938_crop_centers.py) | Compare the original and capped crop centres for signal event 938. |
| [compare_problematic_signal_crops.py](compare_problematic_signal_crops.py) | Compare old/new ROOT crops for the previously problematic signal events. |
| [diagnose_event6008_background.py](diagnose_event6008_background.py) | Show why the old free fit for event 6008 followed background activity. |
| [counterfactual_event6008_new_pilots.py](counterfactual_event6008_new_pilots.py) | Remove eFlag==1 tracks from event 6008 and re-evaluate both new pilots. |
| [evaluate_problematic_new_pilots.py](evaluate_problematic_new_pilots.py) | Evaluate the previously problematic signal events with the new centering pilots. |
| [inspect_truth_signal_tracks.py](inspect_truth_signal_tracks.py) | Inspect signal-only FEDRA segments for selected nue events. |
| [diagnose_background_candidate_z_order.py](diagnose_background_candidate_z_order.py) | Compare z-order sensitivity of scan candidates and matched signal controls. |

These scripts document experiment-specific comparisons and historical debugging. Some require the original CERNBox/EOS inputs or temporary scan-split files. They are preserved for reproducibility, not all intended as current production entry points.
