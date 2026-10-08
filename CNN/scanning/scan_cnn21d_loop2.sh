#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1

brick=$1
volume_dir="/Users/fabioali/cernbox/CNN/data_scan_volumes/${brick}"
prediction_root="/Users/fabioali/cernbox/CNN/data_scan_predictions/${brick}"

for normalization in mu_residual; do
  config="training/configs/cnn21d_signal_p24_l5_fitlt5_${normalization}.yaml"
  run_dir="runs/cnn21d_signal_p24_l5/fitlt5_${normalization}/full"
  output_dir="${prediction_root}/fitlt5_${normalization}"
  mkdir -p "$output_dir"

  for volume in "$volume_dir"/data_scan_*.h5; do
    stem=$(basename "$volume" .h5)

      python3 -m scanning.scan_cnn21d_volumes \
      --input "$volume" \
      --config "$config" \
      --checkpoint "$run_dir/best_classification_model.pt" \
      --regression-checkpoint "$run_dir/best_regression_model.pt" \
      --regression-min-score 0.90 \
      --output "$output_dir/${stem}.predictions.h5" \
      --crop-size 20 \
      --stride 10 \
      --threshold 0.90 \
      --batch-size 64
  done

#    python3 -m candidates.export_background_scan_candidates \
#    --input-glob "$output_dir/*.predictions.h5" \
#    --output-csv "$output_dir/candidates_t090.csv" \
#    --output-json "$output_dir/summary_t090.json" \
#    --thresholds 0.90 0.95 0.98 0.99
done
