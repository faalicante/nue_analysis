# ROOT to HDF5 conversion

## Inspected ROOT schema

The repository contains one ignored ROOT artifact, `event25_samples.root`.  It was
inspected directly rather than inferred from the current C++ source.  On 2026-08-19
its schema was:

- TTree `samples`, title `Raw source regions`, 1 entry;
- `counts`: `int32_t[57][32][32]`, stored as `[z,y,x]`;
- `background_mu`: scalar `float`;
- `presence`: scalar `int32_t`;
- `slope_x`, `slope_y`: scalar `float` in mrad;
- `signal_event_id`, `cell_id`, `tag_cell_id`: scalar `int32_t`.

The entry is a positive sample (`presence=1`, `signal_event_id=25`) with raw-count
range 0--102, no negative counts, and `background_mu=3.7664256`.  The source has
no `sample_type` or negative-type branch.

The current `data_preparation/crop_bw.C` producer and the reference artifact both
use `57 x 32 x 32` crops. The converter requires exactly that shape and exits
with a dimension-mismatch error for other dimensions. It never crops, pads,
reshapes, or resamples a mismatched volume.

The production code uses 50 micrometre xy bins and a 1350 micrometre z step.  The
stored slopes are mrad, so the HDF5 displacement is
`slope_xy = (slope_x, slope_y) * 1350 / (1000 * 50)`, in xy bins per z bin.

## Install and inspect

```sh
python3 -m venv /tmp/crop-bw-venv
/tmp/crop-bw-venv/bin/python -m pip install -r requirements-root-to-hdf5.txt
/tmp/crop-bw-venv/bin/python -m data_preparation.root_to_hdf5 inspect event25_samples.root --json
```

Tree and branch roles are auto-resolved only when unambiguous.  Every role also
has an explicit override such as `--tree`, `--regions-branch`, `--tile-id-branch`,
or `--sample-type-branch`.

## Reduced conversion first

Use `--max-events` for the mandatory reduced validation run:

```sh
/tmp/crop-bw-venv/bin/python -m data_preparation.root_to_hdf5 convert event25_samples.root \
  --output sample_check.h5 --max-events 100
```

The command writes `sample_check.h5`, `sample_check.manifest.csv`,
`sample_check.statistics.json`, `sample_check.schema.json`, and QA projections in
`sample_check_qa/`.  Review those artifacts before running without `--max-events`.

For a complete conversion after review:

```sh
/tmp/crop-bw-venv/bin/python -m data_preparation.root_to_hdf5 convert /path/to/root/files \
  --output dataset.h5 --chunk-size 16 --compression gzip --compression-level 4
```

By default the converter excludes signal entries for which any z slice has a
nonzero support span smaller than `20x20`. This detects a crop lost at a map
edge while ignoring isolated zero-count bins expected from Poisson statistics.
The statistics JSON records every excluded `signal_event_id` and its minimum
support. Use `--min-signal-support 0` only when this filter must be disabled.
Hard and Poisson entries from the known-bad `cell_id=58` are also excluded by
default. Additional background tiles can be excluded by repeating
`--exclude-background-tile TILE_ID`.

Inputs may be ROOT files or directories (searched recursively).  Split-looking
directories or delimited file names (`train`, `validation`/`val`, `test`) preserve
an existing assignment.  Otherwise, a stable hash of `signal_event_id` for signal
or `tile_id`/`cell_id` for background assigns 70/15/15 splits.  Conversion stops
if the same group is found in multiple pre-existing splits.

## HDF5 schema

- `regions_raw`: `[N,57,32,32]`, source integer dtype, `[event,z,y,x]`;
- `background_mu`: semanticamente uno scalare per evento, serializzato come
  `[57]` solo se è comune a tutti gli eventi oppure `[N,57]`; le 57 posizioni di
  ogni evento sono copie identiche, non stime diverse per z;
- `presence`: `[N]`, `uint8`, 0 no shower (Poisson) and 1 shower (hard or signal);
- `sample_type`: `[N]`, 0 Poisson, 1 hard, 2 signal; this is the classification label;
- `slope_xy`: `[N,2]`, `float32`, `(sx,sy)` in xy bins per z bin;
- `signal_event_id`, `tile_id`, `crop_id_in_tile`: `[N]`, with `-1` when absent or
  not applicable;
- `negative_type`: `[N]`, `-1` not applicable, 0 Poisson, 1 hard, 2 unknown;
- `split`: `[N]`, with mapping stored in an attribute;
- `split_indices/{train,validation,test}`: row indices;
- `source_file`, `source_entry`: provenance;
- `metadata/*`: additional finite numeric scalar ROOT branches, when present.

Datasets are chunked and compressed.  Counts are copied verbatim; no normalization,
threshold, smoothing, or augmentation is performed.

## Tests

```sh
/tmp/crop-bw-venv/bin/python -m pytest -q visualization/tests/test_root_to_hdf5.py
```

The synthetic tests cover byte-for-byte count preservation, slope conversion,
hard-negative shower-presence preservation, common/per-event background shapes, grouped
splits, leakage detection, compression/chunking, invalid backgrounds, QA figures,
and dimension mismatch failure.
