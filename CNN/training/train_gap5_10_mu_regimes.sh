#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1


set -euo pipefail

if (( $# < 1 || $# > 2 )); then
    echo "Usage: $0 <pilot|intermediate|full> [low|high|both]"
    exit 1
fi

mode=$1
regime=${2:-both}

case "$mode" in
    pilot|intermediate|full) ;;
    *)
        echo "Invalid mode: $mode"
        exit 1
        ;;
esac

train_regime() {
    local selected_regime=$1
    echo "Training gap5_10 ${selected_regime}-mu model (${mode})"
    python3 -m training.train_cnn21d \
        --config "training/configs/cnn21d_signal_p24_l5_gap5_10_mu_${selected_regime}_full.yaml" \
        --mode "$mode"
}

case "$regime" in
    low|high)
        train_regime "$regime"
        ;;
    both)
        train_regime low
        train_regime high
        ;;
    *)
        echo "Invalid regime: $regime"
        exit 1
        ;;
esac
