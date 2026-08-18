# Source-region export

Build:

```sh
make
```

Run on a positive electron-neutrino event:

```sh
./run.exe 1 EVENT --input input.root --output samples.root \
  --x0 X --y0 Y --tx TX --ty TY --p0 PLATE
```

The output file is opened in `RECREATE` mode, so each run overwrites it with a
new `samples` TTree containing one entry per exported source region. Its branches are:

- `counts[57][20][20]` (`Int_t`), filled in `[z,y,x]` order from the raw,
  unsmoothed per-plate maps;
- `background_mu`, the raw MPV background estimate computed by `crop_bw.C` and
  divided by the 57 plates;
- `presence`, `sample_type`, `slope_x`, `slope_y`, `signal_event_id`, `cell_id`, and
  `tag_cell_id`.

`sample_type=0` identifies Poisson background (`presence=0`), `sample_type=1`
identifies hard negatives (`presence=1`, zero slopes), and `sample_type=2`
identifies signal (`presence=1`, true slopes). Both background classes use
`signal_event_id=-1`. Metadata which is not available is stored as `-1`. Crops
crossing a map boundary are zero-padded; observed bins are copied without
normalization, smoothing, thresholding, or augmentation.

Candidate localization and `zScale` use smoothed copies of the per-plate maps
and their smoothed sum. Across all 57 plates, hard-negative `zScale` uses the
central `10x10` bins while Poisson `zScale` uses the complete `20x20` crop.
Only the independent raw maps are copied into `counts`.
The candidate threshold is derived from the smoothed spectrum, while
`background_mu` stores the raw-spectrum MPV divided by 57.

The corresponding input modes are `data=3` for Poisson background, `data=0`
for hard negatives, and `data=1` for signal.

`--input` overrides the hard-coded production path and `--output` defaults to
`samples.root`. Every candidate passing the existing physics selections is
written; candidate finding examines up to `ntag=500` maxima. PNG generation is
disabled by default; `--images` explicitly enables the legacy image output when
needed for diagnostics.

Render both raw `20x20` slices and smoothed central `20x20` crops directly from
the tree:

```sh
python3 tree_to_images.py event25_samples.root --output event25_tree_images
```

The script uses the scalar `background_mu` as the color minimum of every slice
and one common color maximum across all 57 slices for each output version. It
writes the images under `raw/` and `smooth_crop20/`.


The Makefile contains all of the compilation instructions, the .exe is not necessary here but I find it nicer to have
