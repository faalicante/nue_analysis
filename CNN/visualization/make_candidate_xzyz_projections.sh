#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1

brick=$1
model="$2"

python -m visualization.make_candidate_xzyz_projections \
  --volume-glob "/Users/fabioali/cernbox/CNN/data_scan_volumes/b00$brick/*.h5" \
  --prediction-glob "/Users/fabioali/cernbox/CNN/data_scan_predictions/b00$brick/$model/*.predictions.h5" \
  --output-dir "/Users/fabioali/cernbox/CNN/data_xzyz/b00$brick/$model" \
  --sample-name b000$brick \
  --display-crop-size 40
