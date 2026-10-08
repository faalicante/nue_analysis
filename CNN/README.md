# CNN

SND@LHC shower identification: prepare ROOT/HDF5 samples, train a two-headed 2+1D CNN, scan full cells, and inspect candidates for vertex matching.

## Code categories

Each category has a short README describing every script.

| Folder | Purpose |
|---|---|
| [data_preparation/](data_preparation/README.md) | ROOT crop production, HDF5 conversion, dataset selection and split alignment. |
| [training/](training/README.md) | Dataset/DataLoader, network, loss, metrics, training launchers and configurations. |
| [scanning/](scanning/README.md) | Signal/background/data volume production and sliding-window CNN inference. |
| [candidates/](candidates/README.md) | Candidate export, post-selection, manual comparisons and vertex-search seeds. |
| [studies/](studies/README.md) | Threshold, normalization, efficiency, angle and event-specific studies. |
| [visualization/](visualization/README.md) | GIFs, projections, diagnostic plots, automated tests, utilities and presentation sources. |

## Running the code

Run commands from this `CNN` directory. Python scripts use module syntax so imports work across categories:

```sh
python -m pip install -r requirements-cnn.txt
python -m data_preparation.root_to_hdf5 --help
python -m training.train_cnn21d --help
python -m scanning.scan_cnn21d_volumes --help
python -m pytest -q
```

For an explicit pilot training run:

```sh
python -m training.train_cnn21d --config training/configs/cnn21d.yaml --mode pilot
```

Run shell launchers with `bash category/script.sh`; they locate the project root automatically. Build the C++ producer with `make` (requires ROOT), then use `./data_preparation/run.exe`. Existing script options and scientific selections are unchanged. External CERNBox/EOS paths still require access to the original inputs.

## Preserved data and results

- `runs/`: complete experiment results, checkpoints, metrics and graphics, preserved unchanged.
- `candidate_cut_diagnostics/`: preserved diagnostic reference material.
- Root-level HDF5 datasets and their manifests/schema/statistics: kept at their existing relative paths for compatibility with saved experiments.
- `scan_shards_*`, `scan_signal_*`, `hard_cell58_filter_check.*`, `signal_support_check_100.*`, `converter_retest_samples/` and `event25_samples.root`: preserved test/reference products.
- `manual_showers_b000121.txt`: manual annotations; `outputs/`: candidate notes.
- `output/`: final presentations and supporting figures/tables. The H-mu presentation and asset package retained are version 4; both English and Italian project overviews remain.

Generated caches, compiled binaries, temporary chart exports, superseded slide versions, presentation build dependencies/previews and standalone QA image batches were removed during the October 2026 cleanup. Presentation source code and input snapshots are in `visualization/presentations/`. Historical paths recorded inside preserved results are provenance and were not rewritten.
