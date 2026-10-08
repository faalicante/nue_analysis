#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1

# Scan a separate hard-background ROOT production with the gap5_10 model.

set -euo pipefail

if (( $# != 2 )); then
    echo "Usage: $0 <tag> <root-template>"
    echo "Example template: '/eos/.../cell_{x}0_{y}0/b000021/b000021.0.{x}.{y}.trk.root'"
    exit 1
fi

tag=$1
root_template=$2
model="gap5_10"
volume_dir="/Users/fabioali/cernbox/CNN/high_mu_hard_scan_volumes/$tag"
prediction_dir="/Users/fabioali/cernbox/CNN/high_mu_hard_scan_predictions/$tag/$model"
summary_dir="/Users/fabioali/cernbox/CNN/high_mu_hard_scan_summaries/$tag/$model"

mkdir -p "$prediction_dir" "$summary_dir"

python3 -m scanning.produce_background_scan_shards \
    --root-template "$root_template" \
    --output-dir "$volume_dir" \
    --shard-size 25 \
    --cell-xmin-um 5729 \
    --cell-ymin-um 196206 \
    --skip-errors

for input_path in "$volume_dir"/background_scan_*.h5; do
    [[ -e "$input_path" ]] || continue
    stem=$(basename "$input_path" .h5)
    python3 -m scanning.scan_cnn21d_volumes \
        --input "$input_path" \
        --config "training/configs/cnn21d_signal_p24_l5_${model}_full.yaml" \
        --checkpoint "runs/cnn21d_signal_p24_l5/${model}/full/best_classification_model.pt" \
        --regression-checkpoint "runs/cnn21d_signal_p24_l5/${model}/full/best_regression_model.pt" \
        --output "$prediction_dir/${stem}.predictions.h5" \
        --crop-size 20 \
        --stride 10 \
        --threshold 0.90 \
        --batch-size 64
done

python3 -m candidates.export_background_scan_candidates \
    --input-glob "$prediction_dir/background_scan_*.predictions.h5" \
    --output-csv "$summary_dir/candidates_t090.csv" \
    --output-json "$summary_dir/summary.json" \
    --thresholds 0.80 0.85 0.90 0.95 0.97 0.99
