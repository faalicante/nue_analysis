#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1


set -euo pipefail

if (( $# != 1 )); then
    echo "Usage: $0 <pilot|intermediate|full>"
    exit 1
fi

mode=$1
case "$mode" in
    pilot|intermediate|full) ;;
    *)
        echo "Invalid mode: $mode"
        exit 1
        ;;
esac

python3 -m training.train_cnn21d \
    --config training/configs/cnn21d_signal_p24_l5_gap5_10_mu_aug_global_full.yaml \
    --mode "$mode"

