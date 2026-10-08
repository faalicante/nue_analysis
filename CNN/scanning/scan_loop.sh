#!/usr/bin/env bash

# Resolve project-relative data/configuration paths from any working directory.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1

bricks=(411 222 224)
#gaps=(10)

for brick in "${bricks[@]}"; do
#	for gap in "${gaps[@]}"; do
		echo "scanning/scan_cnn21d_loop.sh $brick 10"
		bash scanning/scan_cnn21d_loop.sh $brick 10
#	done
done

