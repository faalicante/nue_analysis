#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1

bricks=(1021 1022)
models=(gap5_10_mu_high)

for brick in "${bricks[@]}"; do
	for model in "${models[@]}"; do
		bash visualization/make_candidate_xzyz_projections.sh $brick $model
		bash visualization/make_background_candidate_gifs.sh $brick $model
	done
done

