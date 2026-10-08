# Candidate selection and vertex matching

Run Python tools from the `CNN` project root as `python -m candidates.<script_name>` (omit `.py`). Run shell launchers with `bash candidates/<script_name>.sh`.

| Code | Purpose |
|---|---|
| [export_background_scan_candidates.py](export_background_scan_candidates.py) | Aggregate truth-free background scan shards and export manual-check candidates. |
| [export_exclusive_scan_candidates.py](export_exclusive_scan_candidates.py) | Export target scan candidates that do not overlap a reference model. |
| [postfilter_shower_candidates.py](postfilter_shower_candidates.py) | Apply the Poisson-coherent shower-size cut to scan candidates. |
| [export_vertex_search_seeds.py](export_vertex_search_seeds.py) | Export inclined-shower candidates as backward tracks for vertex matching. |
| [validate_shower_start_mc.py](validate_shower_start_mc.py) | Validate the diagnostic shower-start estimator against signal MC metadata. |
| [analyze_manual_data_showers.py](analyze_manual_data_showers.py) | Compare manually identified data showers with every CNN scan window. |

Manual annotations remain at the project root in `manual_showers_b000121.txt`. Candidate tables and manifests should be kept alongside the corresponding predictions.
