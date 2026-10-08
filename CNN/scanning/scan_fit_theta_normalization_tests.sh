#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1


set -euo pipefail

if (( $# < 1 || $# > 2 )); then
    echo "Usage: $0 <pilot|intermediate|full> [signal|background|both]"
    exit 1
fi

mode=$1
sample=${2:-both}
case "$mode" in pilot|intermediate|full) ;; *) echo "Invalid mode: $mode"; exit 1 ;; esac
case "$sample" in signal|background|both) ;; *) echo "Invalid sample: $sample"; exit 1 ;; esac

python_bin=${PYTHON_BIN:-python3}
signal_volume=${SIGNAL_VOLUME:-/private/tmp/crop_bw_scan_splits_l5/signal_scan_test.h5}
background_dir=${BACKGROUND_DIR:-/Users/fabioali/cernbox/CNN/background_scan_volumes}
threshold=0.90

scan_one() {
    local input=$1
    local config=$2
    local classifier=$3
    local regression=$4
    local output=$5
    local report="${output%.h5}.report.json"

    if [[ -s "$output" && -s "$report" ]]; then
        echo "Already complete: $output"
        return
    fi
    mkdir -p "$(dirname "$output")"
    "$python_bin" -m scanning.scan_cnn21d_volumes \
        --input "$input" \
        --config "$config" \
        --checkpoint "$classifier" \
        --regression-checkpoint "$regression" \
        --regression-min-score "$threshold" \
        --output "$output" \
        --crop-size 20 \
        --stride 10 \
        --threshold "$threshold" \
        --batch-size 64 \
        --truth-propagation-limit 40
}

for normalization in mu_residual mpv_z_sigma mpv_z_residual; do
    config="training/configs/cnn21d_signal_p24_l5_fitlt5_${normalization}.yaml"
    run_dir="runs/cnn21d_signal_p24_l5/fitlt5_${normalization}/${mode}"
    classifier="$run_dir/best_classification_model.pt"
    regression="$run_dir/best_regression_model.pt"
    scan_dir="$run_dir/scan_t090"

    if [[ ! -s "$classifier" || ! -s "$regression" ]]; then
        echo "Missing checkpoints for ${normalization} (${mode})" >&2
        exit 1
    fi

    if [[ "$sample" == signal || "$sample" == both ]]; then
        signal_prediction="$scan_dir/signal_scan_test_t090.h5"
        scan_one "$signal_volume" "$config" "$classifier" "$regression" "$signal_prediction"
        "$python_bin" -m studies.summarize_scan_volume_predictions \
            --volumes "$signal_volume" \
            --predictions "$signal_prediction" \
            --output "$scan_dir/signal_scan_summary.json" \
            --thresholds 0.90 0.95 0.98 0.99 \
            --truth-propagation-limit 40
    fi

    if [[ "$sample" == background || "$sample" == both ]]; then
        background_output="$scan_dir/background_scan_t090"
        mkdir -p "$background_output"
        while IFS= read -r input; do
            name=$(basename "$input")
            stem=${name%.h5}
            scan_one "$input" "$config" "$classifier" "$regression" "$background_output/${stem}_t090.h5"
        done < <(find "$background_dir" -maxdepth 1 -type f -name 'background_scan_*.h5' | sort)

        "$python_bin" -m candidates.export_background_scan_candidates \
            --input-glob "$background_output/background_scan_*_t090.h5" \
            --output-csv "$background_output/candidates_t090.csv" \
            --output-json "$background_output/summary.json" \
            --thresholds 0.90 0.95 0.98 0.99

        if [[ "$sample" == both ]]; then
            "$python_bin" -m studies.evaluate_scan_score_size_operating_points \
                --background-volume-glob "$background_dir/background_scan_*.h5" \
                --background-prediction-glob "$background_output/background_scan_*_t090.h5" \
                --signal-volumes "$signal_volume" \
                --signal-predictions "$scan_dir/signal_scan_test_t090.h5" \
                --thresholds 0.90 0.95 0.98 0.99 \
                --theta-min 10 \
                --min-component-voxels 8 \
                --truth-propagation-limit 40 \
                --output-dir "$scan_dir/score_size_operating_points"
        fi
    fi
done
