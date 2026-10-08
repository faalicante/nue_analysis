#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1


set -euo pipefail

python_bin=${PYTHON_BIN:-python3}
signal_volume=${SIGNAL_VOLUME:-/private/tmp/crop_bw_scan_splits_l5/signal_scan_test.h5}
background_root=${BACKGROUND_ROOT:-/Users/fabioali/cernbox/CNN/background_scan_volumes}
threshold=0.90
seed=20260910

run_scan() {
    local input_path=$1
    local config_path=$2
    local classification_checkpoint=$3
    local regression_checkpoint=$4
    local target_mu=$5
    local output_path=$6
    local partial_path="${output_path%.h5}.partial.h5"
    local report_path="${output_path%.h5}.report.json"
    local partial_report_path="${partial_path%.h5}.report.json"

    if [[ -s "$output_path" && -s "$report_path" ]]; then
        echo "Already complete: $output_path"
        return
    fi
    mkdir -p "$(dirname "$output_path")"
    "$python_bin" -m scanning.scan_cnn21d_volumes \
        --input "$input_path" \
        --config "$config_path" \
        --checkpoint "$classification_checkpoint" \
        --regression-checkpoint "$regression_checkpoint" \
        --output "$partial_path" \
        --crop-size 20 \
        --stride 10 \
        --threshold "$threshold" \
        --batch-size 128 \
        --regression-min-score 0.50 \
        --truth-propagation-limit 40 \
        --poisson-background-target-mu "$target_mu" \
        --poisson-background-seed "$seed"
    mv "$partial_path" "$output_path"
    mv "$partial_report_path" "$report_path"
}

summarize_signal() {
    local prediction_path=$1
    local output_path=$2
    "$python_bin" -m studies.summarize_scan_volume_predictions \
        --volumes "$signal_volume" \
        --predictions "$prediction_path" \
        --output "$output_path" \
        --thresholds 0.50 0.90 0.95 0.98 0.99 \
        --truth-propagation-limit 40
}

scan_background_set() {
    local input_dir=$1
    local output_dir=$2
    local config_path=$3
    local classification_checkpoint=$4
    local regression_checkpoint=$5
    local target_mu=$6

    local input_path
    local file_name
    local stem
    while IFS= read -r input_path; do
        file_name=$(basename "$input_path")
        if [[ ! "$file_name" =~ ^background_scan_[0-9]{4}_cell[0-9]{3}-[0-9]{3}\.h5$ ]]; then
            continue
        fi
        stem=${file_name%.h5}
        run_scan \
            "$input_path" "$config_path" \
            "$classification_checkpoint" "$regression_checkpoint" \
            "$target_mu" "$output_dir/${stem}.predictions.h5"
    done < <(find "$input_dir" -maxdepth 1 -type f -name 'background_scan_*.h5' | sort)

    "$python_bin" -m candidates.export_background_scan_candidates \
        --input-glob "$output_dir/background_scan_*.predictions.h5" \
        --output-csv "$output_dir/candidates_t090.csv" \
        --output-json "$output_dir/summary.json" \
        --thresholds 0.50 0.90 0.95 0.98 0.99
}

evaluate_regime() {
    local regime=$1
    local target_mu=$2
    local run_dir="runs/cnn21d_signal_p24_l5/gap5_10_mu_${regime}/full"
    local config_path="training/configs/cnn21d_signal_p24_l5_gap5_10_mu_${regime}_full.yaml"
    local classification_checkpoint="$run_dir/best_classification_model.pt"
    local regression_checkpoint="$run_dir/best_regression_model.pt"
    local evaluation_dir="$run_dir/scan_evaluation"
    local mu_tag=${target_mu/./p}
    local signal_prediction="$evaluation_dir/signal_test_mu${mu_tag}_t090.predictions.h5"

    echo "Evaluating ${regime}-mu network at mu=${target_mu}"
    run_scan \
        "$signal_volume" "$config_path" \
        "$classification_checkpoint" "$regression_checkpoint" \
        "$target_mu" "$signal_prediction"
    summarize_signal "$signal_prediction" "$evaluation_dir/signal_test_mu${mu_tag}_summary.json"

    scan_background_set \
        "$background_root" "$evaluation_dir/background_hard_mu${mu_tag}_t090" \
        "$config_path" "$classification_checkpoint" "$regression_checkpoint" "$target_mu"
    scan_background_set \
        "$background_root/b000024" "$evaluation_dir/background_b000024_mu${mu_tag}_t090" \
        "$config_path" "$classification_checkpoint" "$regression_checkpoint" "$target_mu"
}

evaluate_regime low 5.0
evaluate_regime high 7.5
