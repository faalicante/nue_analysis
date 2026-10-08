#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1


if (( $# !=2 )); then
    echo "Usage: $0 <brick> <gap>"
    #exit 1
fi

brick=$1
gap=$2

INPUT_DIR="/Users/fabioali/cernbox/CNN/data_scan_volumes/b000${brick}"
OUTPUT_DIR="/Users/fabioali/cernbox/CNN/data_scan_predictions/b000${brick}"
MODEL="gap5_${gap}_mu_high"
mkdir -p "$OUTPUT_DIR/$MODEL"

find "$INPUT_DIR" -maxdepth 1 -type f \
    -name "data_scan_b000${brick}_????_*.h5" \
    | sed -E 's/.*_([0-9]{4})_.*\.h5/\1/' \
    | sort -n -u \
    | while IFS= read -r volume
do
    echo "========================================"
    echo "Processing volume $volume"
    echo "========================================"

    python3 -m scanning.scan_cnn21d_volumes \
        --input "$INPUT_DIR/data_scan_b000${brick}_${volume}_"*.h5 \
        --config "training/configs/cnn21d_signal_p24_l5_${MODEL}_full.yaml" \
        --checkpoint "runs/cnn21d_signal_p24_l5/${MODEL}/full/best_classification_model.pt" \
        --regression-checkpoint "runs/cnn21d_signal_p24_l5/${MODEL}/full/best_regression_model.pt" \
        --output "$OUTPUT_DIR/$MODEL/data_scan_b000${brick}_${volume}.predictions.h5" \
        --crop-size 20 \
        --stride 10 \
        --threshold 0.90 \
        --batch-size 64

    if (( $? != 0 )); then
        echo "ERROR processing volume $volume"
        #exit 1
    fi
done

echo
echo "All volumes completed."
