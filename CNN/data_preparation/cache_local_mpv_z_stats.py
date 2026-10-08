#!/usr/bin/env python3
"""Cache per-layer local MPV and lower-side width for a CNN training HDF5."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

from training.cnn_dataset import local_mpv_sigma_per_z


def create_or_resume(
    metadata: h5py.Group, count: int, chunk_size: int, overwrite: bool
) -> tuple[h5py.Dataset, h5py.Dataset, h5py.Dataset]:
    names = ("local_mpv_z", "local_sigma_z", "local_background_stats_complete")
    existing = [name in metadata for name in names]
    if any(existing) and not all(existing):
        raise ValueError("incomplete local-stat cache schema; remove all three datasets before retrying")
    if all(existing) and overwrite:
        for name in names:
            del metadata[name]
        existing = [False, False, False]
    if not any(existing):
        mpv = metadata.create_dataset(
            "local_mpv_z", shape=(count, 57), dtype=np.float32,
            chunks=(min(chunk_size, count), 57), compression="gzip", compression_opts=4, shuffle=True,
        )
        sigma = metadata.create_dataset(
            "local_sigma_z", shape=(count, 57), dtype=np.float32,
            chunks=(min(chunk_size, count), 57), compression="gzip", compression_opts=4, shuffle=True,
        )
        complete = metadata.create_dataset(
            "local_background_stats_complete", shape=(count,), dtype=np.bool_,
            chunks=(min(chunk_size, count),), compression="gzip", compression_opts=4, shuffle=True,
        )
        mpv.attrs["definition"] = "modal raw voxel count in each 32x32 z layer"
        sigma.attrs["definition"] = "lower-side RMS around the local MPV in each 32x32 z layer"
        return mpv, sigma, complete
    mpv, sigma, complete = (metadata[name] for name in names)
    if mpv.shape != (count, 57) or sigma.shape != (count, 57) or complete.shape != (count,):
        raise ValueError("cached local-stat shapes do not match regions_raw")
    return mpv, sigma, complete


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=16)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.chunk_size <= 0:
        raise ValueError("chunk-size must be positive")

    path = args.hdf5.expanduser().resolve()
    with h5py.File(path, "r+") as hdf5:
        raw = hdf5["regions_raw"]
        if raw.shape[1:] != (57, 32, 32):
            raise ValueError(f"expected regions_raw [N,57,32,32], found {raw.shape}")
        metadata = hdf5.require_group("metadata")
        mpv_ds, sigma_ds, complete_ds = create_or_resume(
            metadata, len(raw), args.chunk_size, args.overwrite
        )
        complete = np.asarray(complete_ds[:], dtype=bool)
        pending = np.flatnonzero(~complete)
        print(f"cached={int(complete.sum())} pending={len(pending)} total={len(raw)}", flush=True)
        for first in range(0, len(pending), args.chunk_size):
            indices = pending[first : first + args.chunk_size]
            values = np.asarray(raw[indices])
            mpv = np.empty((len(indices), 57), dtype=np.float32)
            sigma = np.empty((len(indices), 57), dtype=np.float32)
            for row, volume in enumerate(values):
                mpv[row], sigma[row] = local_mpv_sigma_per_z(volume)
            mpv_ds[indices] = mpv
            sigma_ds[indices] = sigma
            complete_ds[indices] = True
            hdf5.flush()
            done = first + len(indices)
            print(f"cached {done}/{len(pending)} pending rows", flush=True)
        hdf5.attrs["local_background_stats"] = "per-z MPV and lower-side RMS on raw 32x32 crop"
        hdf5.attrs["local_background_stats_complete"] = bool(np.all(complete_ds[:]))
        print(f"complete={bool(np.all(complete_ds[:]))}", flush=True)


if __name__ == "__main__":
    main()
