#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1

for input_path in /Users/fabioali/cernbox/CNN/signal_scan_volumes/signal_scan_*.h5; do
    file_name="${input_path##*/}"
    stem="${file_name%.h5}"

    python -m scanning.scan_cnn21d_volumes \
        --input "$input_path" \
        --config training/configs/cnn21d_signal_p24_l5_dual_checkpoint_full_resume.yaml \
        --checkpoint runs/cnn21d_signal_p24_l5/dual_checkpoint_balanced_resume/full/best_classification_model.pt \
        --regression-checkpoint runs/cnn21d_signal_p24_l5/dual_checkpoint_balanced_resume/full/best_regression_model.pt \
        --output "runs/cnn21d_signal_p24_l5/dual_checkpoint_balanced_resume/full/signal_scan_t090/${stem}_t090.h5" \
        --crop-size 20 \
        --stride 10 \
        --threshold 0.90 \
        --batch-size 64
done
