#!/usr/bin/env python3
"""Load one real batch from every HDF5 split and print a JSON summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from training.cnn_dataset import CNNDataConfig, CNNLoaderConfig, create_cnn_dataloaders


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hdf5", type=Path)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--num-workers", type=int, default=2)
    parser.add_argument("--batches", type=int, default=2)
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()

    loaders = create_cnn_dataloaders(
        args.hdf5,
        data_config=CNNDataConfig(max_jitter=5, seed=args.seed),
        loader_config=CNNLoaderConfig(
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            shuffle_train=False,
            seed=args.seed,
        ),
    )
    summary = {}
    for split, loader in loaders.items():
        iterator = iter(loader)
        batches = [next(iterator) for _ in range(args.batches)]
        batch = batches[0]
        summary[split] = {
            "volume": list(batch["volume"].shape),
            "volume_dtype": str(batch["volume"].dtype),
            "sample_type": list(batch["sample_type"].shape),
            "sample_type_dtype": str(batch["sample_type"].dtype),
            "slope_xy": list(batch["slope_xy"].shape),
            "slope_dtype": str(batch["slope_xy"].dtype),
            "offset_min_xy": batch["metadata"]["crop_offset_xy"].amin(dim=0).tolist(),
            "offset_max_xy": batch["metadata"]["crop_offset_xy"].amax(dim=0).tolist(),
            "rotation_k": sorted(
                set(batch["metadata"]["rotation_quarter_turns_ccw"].tolist())
            ),
            "worker_ids": sorted(
                {
                    worker_id
                    for current in batches
                    for worker_id in current["metadata"]["worker_id"].tolist()
                }
            ),
            "worker_pids": sorted(
                {
                    worker_pid
                    for current in batches
                    for worker_pid in current["metadata"]["worker_pid"].tolist()
                }
            ),
        }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
