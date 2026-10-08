#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1

brick=$1
model="$2"

python -m visualization.make_background_candidate_gifs \
  --volume-glob "/Users/fabioali/cernbox/CNN/data_scan_volumes/b00$brick/*.h5" \
  --prediction-glob "/Users/fabioali/cernbox/CNN/data_scan_predictions/b00$brick/$model/*.predictions.h5" \
  --output-dir "/Users/fabioali/cernbox/CNN/data_gifs/b00$brick/$model" \
  --display-crop-size 40 \
  --frame-duration-ms 160
