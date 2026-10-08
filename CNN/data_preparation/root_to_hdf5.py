#!/usr/bin/env python3
"""Convert raw ROOT source-region TTrees to a validated HDF5 dataset.

The converter deliberately does not normalize, threshold, smooth, or augment the
voxel counts.  ROOT volumes are copied in their existing [z, y, x] order.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import h5py
import numpy as np
import uproot


EXPECTED_REGION_SHAPE = (57, 32, 32)
SPLIT_NAMES = ("train", "validation", "test")
SPLIT_TO_CODE = {name: index for index, name in enumerate(SPLIT_NAMES)}
SAMPLE_TYPE_NAMES = {0: "poisson", 1: "hard", 2: "signal"}
NEGATIVE_TYPE_NAMES = {-1: "not_applicable", 0: "poisson", 1: "hard", 2: "unknown"}
DEFAULT_EXCLUDED_BACKGROUND_TILES = (58,)

BRANCH_ALIASES = {
    "regions": ("counts", "regions_raw"),
    "background": ("background_mu",),
    "presence": ("presence",),
    "slope_x": ("slope_x", "sx"),
    "slope_y": ("slope_y", "sy"),
    "signal_event_id": ("signal_event_id", "event_id"),
    "tile_id": ("tile_id", "cell_id"),
    "crop_id": ("crop_id_in_tile", "tag_cell_id"),
    "sample_type": ("sample_type",),
}


class ConversionError(RuntimeError):
    """Raised for a schema or data problem that must stop conversion."""


@dataclass(frozen=True)
class BranchInfo:
    name: str
    typename: str
    interpretation: str
    dtype: str | None
    per_entry_shape: tuple[int, ...] | None


@dataclass
class RootFileSchema:
    path: Path
    tree_name: str
    entries: int
    selected_entries: int
    branches: dict[str, BranchInfo]
    branch_map: dict[str, str | None]
    explicit_split: str | None

    @property
    def regions_dtype(self) -> np.dtype:
        name = self.branch_map["regions"]
        assert name is not None
        dtype = self.branches[name].dtype
        assert dtype is not None
        return np.dtype(dtype)


@dataclass
class Metadata:
    background_mu: np.ndarray
    presence: np.ndarray
    sample_type: np.ndarray
    slope_xy: np.ndarray
    signal_event_id: np.ndarray
    tile_id: np.ndarray
    crop_id_in_tile: np.ndarray
    negative_type: np.ndarray
    source_file: np.ndarray
    source_entry: np.ndarray
    source_split: np.ndarray
    extras: dict[str, np.ndarray]

    def __len__(self) -> int:
        return int(self.presence.shape[0])


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Cannot JSON-encode {type(value).__name__}")


def _expand_inputs(inputs: Sequence[Path]) -> list[Path]:
    expanded: list[Path] = []
    for raw_path in inputs:
        path = raw_path.expanduser()
        if path.is_dir():
            expanded.extend(sorted(candidate for candidate in path.rglob("*.root") if candidate.is_file()))
        elif path.is_file():
            expanded.append(path)
        else:
            raise ConversionError(f"Input does not exist: {path}")

    unique = sorted({path.resolve() for path in expanded}, key=lambda item: str(item))
    if not unique:
        raise ConversionError("No ROOT input files were found")
    return unique


def _infer_existing_split(path: Path) -> str | None:
    # Directory names such as ``train/`` are authoritative.  For a file, accept
    # delimited forms such as ``sample_validation.root``.  Do not treat arbitrary
    # parent names such as ``test_conversion`` as a split marker.
    directory_tokens = {part.lower() for part in path.parent.parts}
    stem_tokens = {token for token in re.split(r"[._-]+", path.stem.lower()) if token}
    tokens = directory_tokens | stem_tokens
    matches: set[str] = set()
    if tokens & {"train", "training"}:
        matches.add("train")
    if tokens & {"val", "valid", "validation"}:
        matches.add("validation")
    if "test" in tokens:
        matches.add("test")
    if len(matches) > 1:
        raise ConversionError(f"Ambiguous existing split in path {path}: {sorted(matches)}")
    return next(iter(matches), None)


def _select_tree(root_file: uproot.ReadOnlyDirectory, requested: str) -> str:
    classnames = root_file.classnames(recursive=False, cycle=False)
    trees = sorted(name for name, classname in classnames.items() if classname == "TTree")
    if requested != "auto":
        if requested not in trees:
            raise ConversionError(
                f"Requested TTree {requested!r} not found; TTrees present: {trees or 'none'}"
            )
        return requested
    if len(trees) != 1:
        raise ConversionError(
            f"Expected exactly one TTree for automatic selection, found {trees or 'none'}; "
            "use --tree"
        )
    return trees[0]


def _sample_branch(branch: Any, entries: int) -> tuple[np.dtype | None, tuple[int, ...] | None]:
    if entries <= 0:
        return None, None
    try:
        values = np.asarray(branch.array(entry_start=0, entry_stop=1, library="np"))
    except Exception:
        return None, None
    if values.shape[0] != 1:
        return values.dtype, None
    return values.dtype, tuple(int(size) for size in values.shape[1:])


def _resolve_branch(
    branch_names: set[str], requested: str, aliases: Sequence[str], role: str, required: bool
) -> str | None:
    if requested != "auto":
        if requested not in branch_names:
            raise ConversionError(f"Branch {requested!r} requested for {role} is absent")
        return requested
    matches = [name for name in aliases if name in branch_names]
    if len(matches) > 1:
        raise ConversionError(f"Multiple branches match {role}: {matches}; select one explicitly")
    if matches:
        return matches[0]
    if required:
        raise ConversionError(f"No branch found for {role}; tried {list(aliases)}")
    return None


def _branch_requests(args: argparse.Namespace) -> dict[str, str]:
    return {
        "regions": args.regions_branch,
        "background": args.background_branch,
        "presence": args.presence_branch,
        "slope_x": args.slope_x_branch,
        "slope_y": args.slope_y_branch,
        "signal_event_id": args.signal_event_id_branch,
        "tile_id": args.tile_id_branch,
        "crop_id": args.crop_id_branch,
        "sample_type": args.sample_type_branch,
    }


def inspect_root_file(
    path: Path, args: argparse.Namespace, selected_entries: int | None = None
) -> RootFileSchema:
    try:
        root_file = uproot.open(path)
    except Exception as error:
        raise ConversionError(f"Cannot open ROOT file {path}: {error}") from error

    with root_file:
        tree_name = _select_tree(root_file, args.tree)
        tree = root_file[tree_name]
        entries = int(tree.num_entries)
        if entries == 0:
            raise ConversionError(f"TTree {tree_name!r} in {path} has no entries")

        branches: dict[str, BranchInfo] = {}
        for name, branch in tree.items():
            dtype, shape = _sample_branch(branch, entries)
            branches[name] = BranchInfo(
                name=name,
                typename=str(getattr(branch, "typename", "")),
                interpretation=repr(getattr(branch, "interpretation", "")),
                dtype=None if dtype is None else dtype.str,
                per_entry_shape=shape,
            )

        branch_names = set(branches)
        requests = _branch_requests(args)
        required_roles = {"regions", "background", "presence", "slope_x", "slope_y"}
        branch_map = {
            role: _resolve_branch(
                branch_names,
                requests[role],
                BRANCH_ALIASES[role],
                role,
                role in required_roles,
            )
            for role in requests
        }

        regions_name = branch_map["regions"]
        assert regions_name is not None
        regions_info = branches[regions_name]
        if regions_info.per_entry_shape != EXPECTED_REGION_SHAPE:
            raise ConversionError(
                f"Dimension mismatch in {path}:{tree_name}/{regions_name}: found "
                f"{regions_info.per_entry_shape}, required {EXPECTED_REGION_SHAPE} in [z,y,x] order. "
                "The converter will not crop, pad, reshape, or resample silently."
            )
        if regions_info.dtype is None or np.dtype(regions_info.dtype).kind not in "iu":
            raise ConversionError(
                f"Raw counts in {path}:{regions_name} must have an integer dtype, found "
                f"{regions_info.dtype}"
            )

        expected_shapes = {
            "background": {(), (EXPECTED_REGION_SHAPE[0],)},
            "presence": {()},
            "slope_x": {()},
            "slope_y": {()},
            "signal_event_id": {()},
            "tile_id": {()},
            "crop_id": {()},
            "sample_type": {()},
        }
        for role, allowed in expected_shapes.items():
            name = branch_map[role]
            if name is not None and branches[name].per_entry_shape not in allowed:
                raise ConversionError(
                    f"Unexpected shape for {role} branch {path}:{name}: "
                    f"{branches[name].per_entry_shape}; allowed {sorted(allowed, key=len)}"
                )

        numeric_roles = {
            "background": "biuf",
            "presence": "biu",
            "slope_x": "biuf",
            "slope_y": "biuf",
            "signal_event_id": "biu",
            "tile_id": "biu",
            "crop_id": "biu",
            "sample_type": "biu",
        }
        for role, allowed_kinds in numeric_roles.items():
            name = branch_map[role]
            if name is None:
                continue
            dtype_string = branches[name].dtype
            if dtype_string is None or np.dtype(dtype_string).kind not in allowed_kinds:
                raise ConversionError(
                    f"Branch {path}:{name} for {role} has unsupported dtype {dtype_string}"
                )

        selected = entries if selected_entries is None else min(entries, selected_entries)
        return RootFileSchema(
            path=path,
            tree_name=tree_name,
            entries=entries,
            selected_entries=selected,
            branches=branches,
            branch_map=branch_map,
            explicit_split=_infer_existing_split(path),
        )


def discover_inputs(paths: Sequence[Path], args: argparse.Namespace) -> list[RootFileSchema]:
    files = _expand_inputs(paths)
    remaining = args.max_events if getattr(args, "max_events", None) is not None else None
    schemas: list[RootFileSchema] = []
    for path in files:
        if remaining is not None and remaining <= 0:
            break
        schema = inspect_root_file(path, args, remaining)
        if schemas and schema.regions_dtype != schemas[0].regions_dtype:
            raise ConversionError(
                f"Inconsistent raw-count dtype: {schemas[0].path} has {schemas[0].regions_dtype}, "
                f"but {schema.path} has {schema.regions_dtype}"
            )
        schemas.append(schema)
        if remaining is not None:
            remaining -= schema.selected_entries
    if not schemas or sum(schema.selected_entries for schema in schemas) == 0:
        raise ConversionError("No entries selected for conversion")
    return schemas


def _read_branch(tree: Any, name: str, stop: int) -> np.ndarray:
    values = np.asarray(tree[name].array(entry_start=0, entry_stop=stop, library="np"))
    if values.shape[0] != stop:
        raise ConversionError(f"Branch {name} yielded {values.shape[0]} entries, expected {stop}")
    return values


def _canonical_ids(values: np.ndarray | None, count: int) -> np.ndarray:
    if values is None:
        return np.full(count, -1, dtype=np.int64)
    result = np.asarray(values, dtype=np.int64).reshape(count)
    result[result < 0] = -1
    return result


def _extra_branch_specs(schemas: Sequence[RootFileSchema]) -> dict[str, np.dtype]:
    specs: dict[str, np.dtype] = {}
    for schema in schemas:
        used = {name for name in schema.branch_map.values() if name is not None}
        for name, info in schema.branches.items():
            if name in used or info.per_entry_shape != () or info.dtype is None:
                continue
            dtype = np.dtype(info.dtype)
            # Boolean optionals cannot represent the -1 missing sentinel without
            # changing meaning, so only preserve scalar integer/float metadata.
            if dtype.kind not in "iuf":
                continue
            if name in specs and specs[name] != dtype:
                raise ConversionError(
                    f"Extra scalar metadata branch {name!r} has inconsistent dtypes: "
                    f"{specs[name]} and {dtype}"
                )
            specs[name] = dtype
    return specs


def read_metadata(
    schemas: Sequence[RootFileSchema], args: argparse.Namespace
) -> tuple[Metadata, list[str]]:
    extra_specs = _extra_branch_specs(schemas)
    columns: dict[str, list[np.ndarray]] = {
        "background_mu": [],
        "presence": [],
        "sample_type": [],
        "slope_xy": [],
        "signal_event_id": [],
        "tile_id": [],
        "crop_id_in_tile": [],
        "negative_type": [],
        "source_file": [],
        "source_entry": [],
        "source_split": [],
    }
    extra_columns: dict[str, list[np.ndarray]] = {name: [] for name in extra_specs}
    warnings: list[str] = []

    slope_factor = args.z_bin_size_um / (args.slope_unit_divisor * args.xy_bin_size_um)
    for schema in schemas:
        count = schema.selected_entries
        with uproot.open(schema.path) as root_file:
            tree = root_file[schema.tree_name]
            mapped = schema.branch_map

            background_name = mapped["background"]
            presence_name = mapped["presence"]
            slope_x_name = mapped["slope_x"]
            slope_y_name = mapped["slope_y"]
            assert background_name and presence_name and slope_x_name and slope_y_name

            background = _read_branch(tree, background_name, count).astype(np.float64, copy=False)
            if background.ndim == 1:
                background = np.repeat(background[:, np.newaxis], EXPECTED_REGION_SHAPE[0], axis=1)
            else:
                background = background.reshape(count, EXPECTED_REGION_SHAPE[0])

            source_presence = _read_branch(tree, presence_name, count).astype(np.int64).reshape(count)
            invalid_source_labels = sorted(set(np.unique(source_presence).tolist()) - {0, 1})
            if invalid_source_labels:
                raise ConversionError(
                    f"Invalid labels in {schema.path}:{presence_name}: {invalid_source_labels}; expected 0 or 1"
                )

            slope_x = _read_branch(tree, slope_x_name, count).astype(np.float64).reshape(count)
            slope_y = _read_branch(tree, slope_y_name, count).astype(np.float64).reshape(count)
            slope_xy = np.column_stack((slope_x, slope_y)) * slope_factor

            signal_name = mapped["signal_event_id"]
            tile_name = mapped["tile_id"]
            crop_name = mapped["crop_id"]
            signal_ids = _canonical_ids(None if signal_name is None else _read_branch(tree, signal_name, count), count)
            tile_ids = _canonical_ids(None if tile_name is None else _read_branch(tree, tile_name, count), count)
            crop_ids = _canonical_ids(None if crop_name is None else _read_branch(tree, crop_name, count), count)

            sample_type_name = mapped["sample_type"]
            if sample_type_name is not None:
                sample_type = _read_branch(tree, sample_type_name, count).astype(np.int64).reshape(count)
                invalid_types = sorted(set(np.unique(sample_type).tolist()) - {0, 1, 2})
                if invalid_types:
                    raise ConversionError(
                        f"Invalid sample_type values in {schema.path}:{sample_type_name}: {invalid_types}"
                    )
                expected_presence = (sample_type != 0).astype(np.int64)
                incompatible = expected_presence != source_presence
                if np.any(incompatible):
                    bad_rows = np.flatnonzero(incompatible)[:10].tolist()
                    raise ConversionError(
                        f"Inconsistent presence/sample_type labels in {schema.path} at rows {bad_rows}"
                    )
            else:
                # Legacy files without sample_type can still be resolved from the
                # shower-presence flag and signal_event_id convention.
                sample_type = np.where(
                    source_presence == 0,
                    0,
                    np.where(signal_ids >= 0, 2, 1),
                ).astype(np.int64)
            presence = source_presence.astype(np.uint8)
            negative_type = np.where(sample_type == 2, -1, sample_type).astype(np.int8)

            if not np.isfinite(background).all():
                raise ConversionError(f"NaN or infinite background_mu found in {schema.path}")
            if np.any(background <= 0):
                minimum = float(np.min(background))
                raise ConversionError(f"background_mu <= 0 found in {schema.path} (minimum {minimum})")
            if not np.isfinite(slope_xy).all():
                raise ConversionError(f"NaN or infinite slope found in {schema.path}")

            signal = sample_type == 2
            background_class = ~signal
            if np.any(signal_ids[signal] < 0):
                raise ConversionError(
                    f"Signal samples in {schema.path} require a non-negative signal_event_id for grouped splitting"
                )
            if np.any(signal_ids[background_class] >= 0):
                raise ConversionError(f"Background samples in {schema.path} have an applicable signal_event_id")
            if np.any(tile_ids[background_class] < 0):
                raise ConversionError(
                    f"Background samples in {schema.path} require a non-negative tile_id/cell_id "
                    "for grouped splitting"
                )
            if np.any(np.all(slope_xy[signal] == 0, axis=1)):
                raise ConversionError(f"Zero two-component slope found for an inclined signal in {schema.path}")
            if np.any(slope_xy[background_class] != 0):
                raise ConversionError(f"Non-zero slope found for background samples in {schema.path}")

            columns["background_mu"].append(background)
            columns["presence"].append(presence)
            columns["sample_type"].append(sample_type.astype(np.int8))
            columns["slope_xy"].append(slope_xy.astype(np.float32))
            columns["signal_event_id"].append(signal_ids)
            columns["tile_id"].append(tile_ids)
            columns["crop_id_in_tile"].append(crop_ids)
            columns["negative_type"].append(negative_type)
            columns["source_file"].append(np.full(count, str(schema.path), dtype=object))
            columns["source_entry"].append(np.arange(count, dtype=np.int64))
            columns["source_split"].append(np.full(count, schema.explicit_split, dtype=object))

            used = {name for name in mapped.values() if name is not None}
            for name, dtype in extra_specs.items():
                if name in schema.branches and name not in used:
                    values = _read_branch(tree, name, count).reshape(count).astype(dtype, copy=False)
                    if dtype.kind == "f" and not np.isfinite(values).all():
                        raise ConversionError(f"NaN or infinite metadata value in {schema.path}:{name}")
                else:
                    # The source does not provide this optional scalar for these
                    # rows.  Keep the HDF5 finite and use the documented missing
                    # metadata sentinel rather than manufacturing NaNs.
                    values = np.full(count, -1, dtype=dtype)
                extra_columns[name].append(values)

    metadata = Metadata(
        background_mu=np.concatenate(columns["background_mu"]),
        presence=np.concatenate(columns["presence"]),
        sample_type=np.concatenate(columns["sample_type"]),
        slope_xy=np.concatenate(columns["slope_xy"]),
        signal_event_id=np.concatenate(columns["signal_event_id"]),
        tile_id=np.concatenate(columns["tile_id"]),
        crop_id_in_tile=np.concatenate(columns["crop_id_in_tile"]),
        negative_type=np.concatenate(columns["negative_type"]),
        source_file=np.concatenate(columns["source_file"]),
        source_entry=np.concatenate(columns["source_entry"]),
        source_split=np.concatenate(columns["source_split"]),
        extras={name: np.concatenate(parts) for name, parts in extra_columns.items()},
    )
    return metadata, warnings


def _subset_metadata(metadata: Metadata, keep: np.ndarray) -> Metadata:
    if keep.shape != (len(metadata),):
        raise ConversionError(
            f"Internal filter-mask shape mismatch: found {keep.shape}, expected {(len(metadata),)}"
        )
    return Metadata(
        background_mu=metadata.background_mu[keep],
        presence=metadata.presence[keep],
        sample_type=metadata.sample_type[keep],
        slope_xy=metadata.slope_xy[keep],
        signal_event_id=metadata.signal_event_id[keep],
        tile_id=metadata.tile_id[keep],
        crop_id_in_tile=metadata.crop_id_in_tile[keep],
        negative_type=metadata.negative_type[keep],
        source_file=metadata.source_file[keep],
        source_entry=metadata.source_entry[keep],
        source_split=metadata.source_split[keep],
        extras={name: values[keep] for name, values in metadata.extras.items()},
    )


def _crop_support_metrics(regions: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return minimum occupied x/y spans and empty-slice counts per event.

    Isolated zero-count bins are expected from Poisson statistics and do not
    reduce the span.  A completely empty edge band, however, moves the first or
    last occupied row/column inward and is therefore treated as lost support.
    """
    occupied = regions != 0
    occupied_x = np.any(occupied, axis=2)  # event,z,x
    occupied_y = np.any(occupied, axis=3)  # event,z,y

    slice_has_content = np.any(occupied_x, axis=2)
    first_x = np.argmax(occupied_x, axis=2)
    last_x_exclusive = EXPECTED_REGION_SHAPE[2] - np.argmax(occupied_x[..., ::-1], axis=2)
    first_y = np.argmax(occupied_y, axis=2)
    last_y_exclusive = EXPECTED_REGION_SHAPE[1] - np.argmax(occupied_y[..., ::-1], axis=2)
    widths = np.where(slice_has_content, last_x_exclusive - first_x, 0)
    heights = np.where(slice_has_content, last_y_exclusive - first_y, 0)
    return (
        np.min(widths, axis=1),
        np.min(heights, axis=1),
        np.count_nonzero(~slice_has_content, axis=1),
    )


def build_signal_support_filter(
    schemas: Sequence[RootFileSchema], metadata: Metadata, args: argparse.Namespace
) -> tuple[list[np.ndarray], dict[str, Any]]:
    minimum = int(args.min_signal_support)
    keep_masks: list[np.ndarray] = []
    rejected: list[dict[str, int | str]] = []
    global_offset = 0

    for schema in schemas:
        count = schema.selected_entries
        keep = np.ones(count, dtype=bool)
        sample_types = metadata.sample_type[global_offset : global_offset + count]
        event_ids = metadata.signal_event_id[global_offset : global_offset + count]
        regions_name = schema.branch_map["regions"]
        assert regions_name is not None

        if minimum > 0 and np.any(sample_types == 2):
            with uproot.open(schema.path) as root_file:
                tree = root_file[schema.tree_name]
                for start in range(0, count, args.chunk_size):
                    stop = min(start + args.chunk_size, count)
                    values = np.asarray(
                        tree[regions_name].array(entry_start=start, entry_stop=stop, library="np")
                    )
                    expected = (stop - start, *EXPECTED_REGION_SHAPE)
                    if values.shape != expected:
                        raise ConversionError(
                            f"Inconsistent chunk shape in {schema.path}:{regions_name}: "
                            f"found {values.shape}, expected {expected}"
                        )
                    if values.dtype.kind not in "iu":
                        raise ConversionError(
                            f"Non-integer raw counts in {schema.path}:{regions_name}"
                        )
                    if np.any(values < 0):
                        minimum_count = int(np.min(values))
                        raise ConversionError(
                            f"Negative raw counts in {schema.path}:{regions_name} "
                            f"(minimum {minimum_count})"
                        )

                    signal = sample_types[start:stop] == 2
                    if not np.any(signal):
                        continue
                    min_width, min_height, empty_slices = _crop_support_metrics(values)
                    unusable = signal & ((min_width < minimum) | (min_height < minimum))
                    for local in np.flatnonzero(unusable):
                        source_entry = start + int(local)
                        keep[source_entry] = False
                        rejected.append(
                            {
                                "source_file": str(schema.path),
                                "source_entry": source_entry,
                                "signal_event_id": int(event_ids[source_entry]),
                                "minimum_nonzero_width": int(min_width[local]),
                                "minimum_nonzero_height": int(min_height[local]),
                                "empty_slices": int(empty_slices[local]),
                            }
                        )

        keep_masks.append(keep)
        global_offset += count

    if global_offset != len(metadata):
        raise ConversionError(
            f"Internal source-row mismatch: inspected {global_offset}, metadata has {len(metadata)}"
        )
    kept = int(sum(np.count_nonzero(mask) for mask in keep_masks))
    report: dict[str, Any] = {
        "enabled": minimum > 0,
        "minimum_signal_support_xy": minimum,
        "criterion": (
            "for every z slice, the bounding span of nonzero raw bins must be at least "
            f"{minimum}x{minimum}; isolated Poisson zeros are ignored"
            if minimum > 0
            else "disabled"
        ),
        "input_entries": len(metadata),
        "kept_entries": kept,
        "excluded_signal_entries": len(rejected),
        "excluded_signal_event_ids": sorted(
            {int(row["signal_event_id"]) for row in rejected}
        ),
        "excluded": rejected,
    }
    return keep_masks, report


def apply_background_tile_filter(
    schemas: Sequence[RootFileSchema],
    metadata: Metadata,
    keep_masks: Sequence[np.ndarray],
    excluded_tiles: Sequence[int],
) -> dict[str, Any]:
    excluded = sorted({int(tile) for tile in excluded_tiles})
    rejected: list[dict[str, int | str]] = []
    global_offset = 0
    for schema, keep in zip(schemas, keep_masks, strict=True):
        count = schema.selected_entries
        sample_types = metadata.sample_type[global_offset : global_offset + count]
        tile_ids = metadata.tile_id[global_offset : global_offset + count]
        crop_ids = metadata.crop_id_in_tile[global_offset : global_offset + count]
        reject = np.isin(sample_types, (0, 1)) & np.isin(tile_ids, excluded)
        for source_entry in np.flatnonzero(reject):
            keep[int(source_entry)] = False
            rejected.append(
                {
                    "source_file": str(schema.path),
                    "source_entry": int(source_entry),
                    "tile_id": int(tile_ids[source_entry]),
                    "crop_id_in_tile": int(crop_ids[source_entry]),
                    "sample_type": int(sample_types[source_entry]),
                    "sample_type_name": SAMPLE_TYPE_NAMES[int(sample_types[source_entry])],
                }
            )
        global_offset += count

    return {
        "enabled": bool(excluded),
        "excluded_background_tile_ids": excluded,
        "excluded_background_entries": len(rejected),
        "excluded_by_class": {
            name: sum(int(row["sample_type"]) == code for row in rejected)
            for code, name in SAMPLE_TYPE_NAMES.items()
            if code in (0, 1)
        },
        "excluded": rejected,
    }


def _hashed_split(kind: str, group_id: int, seed: int, ratios: Sequence[float]) -> str:
    digest = hashlib.blake2b(f"{seed}:{kind}:{group_id}".encode("utf-8"), digest_size=8).digest()
    value = int.from_bytes(digest, "big") / float(1 << 64)
    if value < ratios[0]:
        return "train"
    if value < ratios[0] + ratios[1]:
        return "validation"
    return "test"


def assign_splits(metadata: Metadata, seed: int, ratios: Sequence[float]) -> np.ndarray:
    if len(ratios) != 3 or any(value < 0 for value in ratios) or not np.isclose(sum(ratios), 1.0):
        raise ConversionError(f"Split ratios must be three non-negative values summing to 1, found {ratios}")

    groups: list[tuple[str, int]] = []
    for index, sample_type in enumerate(metadata.sample_type):
        if sample_type == 2:
            groups.append(("signal_event_id", int(metadata.signal_event_id[index])))
        elif sample_type in (0, 1):
            groups.append(("tile_id", int(metadata.tile_id[index])))
        else:
            raise ConversionError(f"Invalid output sample_type {sample_type} at row {index}")

    explicit_by_group: dict[tuple[str, int], str] = {}
    for group, existing in zip(groups, metadata.source_split, strict=True):
        if existing is None:
            continue
        previous = explicit_by_group.get(group)
        if previous is not None and previous != existing:
            raise ConversionError(
                f"Existing split leakage for {group[0]}={group[1]}: {previous} and {existing}"
            )
        explicit_by_group[group] = str(existing)

    assignment_by_group = {
        group: explicit_by_group.get(group, _hashed_split(group[0], group[1], seed, ratios))
        for group in set(groups)
    }
    split = np.asarray([SPLIT_TO_CODE[assignment_by_group[group]] for group in groups], dtype=np.int8)
    validate_group_leakage(metadata, split)
    return split


def validate_group_leakage(metadata: Metadata, split: np.ndarray) -> None:
    seen: dict[tuple[str, int], int] = {}
    for index, (sample_type, split_code) in enumerate(zip(metadata.sample_type, split, strict=True)):
        if int(split_code) not in range(len(SPLIT_NAMES)):
            raise ConversionError(f"Invalid split code {split_code} at row {index}")
        group = (
            ("signal_event_id", int(metadata.signal_event_id[index]))
            if sample_type == 2
            else ("tile_id", int(metadata.tile_id[index]))
        )
        previous = seen.setdefault(group, int(split_code))
        if previous != int(split_code):
            raise ConversionError(
                f"Group leakage: {group[0]}={group[1]} occurs in "
                f"{SPLIT_NAMES[previous]} and {SPLIT_NAMES[int(split_code)]}"
            )


def compute_statistics(metadata: Metadata, split: np.ndarray) -> dict[str, Any]:
    statistics: dict[str, Any] = {
        "total": len(metadata),
        "class": {
            name: int(np.count_nonzero(metadata.sample_type == code))
            for code, name in SAMPLE_TYPE_NAMES.items()
        },
        "presence": {
            "no_shower": int(np.count_nonzero(metadata.presence == 0)),
            "shower": int(np.count_nonzero(metadata.presence == 1)),
        },
        "split": {},
        "negative_type": {},
    }
    for code, name in enumerate(SPLIT_NAMES):
        mask = split == code
        statistics["split"][name] = {
            "total": int(np.count_nonzero(mask)),
            **{
                class_name: int(np.count_nonzero(mask & (metadata.sample_type == class_code)))
                for class_code, class_name in SAMPLE_TYPE_NAMES.items()
            },
        }
    for code, name in NEGATIVE_TYPE_NAMES.items():
        statistics["negative_type"][name] = int(np.count_nonzero(metadata.negative_type == code))
    return statistics


def _dataset_kwargs(
    chunks: tuple[int, ...], compression: str, compression_level: int
) -> dict[str, Any]:
    result: dict[str, Any] = {"chunks": chunks}
    if compression == "none":
        return result
    result.update(compression=compression, shuffle=True)
    if compression == "gzip":
        result["compression_opts"] = compression_level
    return result


def _safe_extra_name(name: str) -> str:
    return name.replace("/", "__")


def write_hdf5(
    output: Path,
    schemas: Sequence[RootFileSchema],
    metadata: Metadata,
    split: np.ndarray,
    keep_masks: Sequence[np.ndarray],
    support_filter: dict[str, Any],
    background_tile_filter: dict[str, Any],
    args: argparse.Namespace,
) -> bool:
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not args.overwrite:
        raise ConversionError(f"Output exists: {output}; use --overwrite to replace it")

    total = len(metadata)
    chunk_rows = min(args.chunk_size, total)
    counts_dtype = schemas[0].regions_dtype
    common_background = all(
        np.array_equal(metadata.background_mu[0], row) for row in metadata.background_mu[1:]
    )

    handle = tempfile.NamedTemporaryFile(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent, delete=False
    )
    temporary = Path(handle.name)
    handle.close()
    try:
        with h5py.File(temporary, "w") as hdf5:
            hdf5.attrs["schema_version"] = "1.0"
            hdf5.attrs["axis_order"] = "event,z,y,x"
            hdf5.attrs["regions_are_raw"] = True
            hdf5.attrs["normalization"] = "none"
            hdf5.attrs["smoothing"] = "none"
            hdf5.attrs["augmentation"] = "none"
            hdf5.attrs["xy_bin_size_um"] = args.xy_bin_size_um
            hdf5.attrs["z_bin_size_um"] = args.z_bin_size_um
            hdf5.attrs["source_slope_unit_divisor"] = args.slope_unit_divisor
            hdf5.attrs["slope_xy_unit"] = "xy_bins_per_z_bin"
            hdf5.attrs["split_names_json"] = json.dumps(SPLIT_NAMES)
            hdf5.attrs["negative_type_names_json"] = json.dumps(NEGATIVE_TYPE_NAMES)
            hdf5.attrs["minimum_signal_support_xy"] = int(
                support_filter["minimum_signal_support_xy"]
            )
            hdf5.attrs["excluded_unusable_signal_entries"] = int(
                support_filter["excluded_signal_entries"]
            )
            hdf5.attrs["excluded_background_tile_ids_json"] = json.dumps(
                background_tile_filter["excluded_background_tile_ids"]
            )
            hdf5.attrs["excluded_background_entries"] = int(
                background_tile_filter["excluded_background_entries"]
            )

            regions = hdf5.create_dataset(
                "regions_raw",
                shape=(total, *EXPECTED_REGION_SHAPE),
                dtype=counts_dtype,
                **_dataset_kwargs(
                    (chunk_rows, *EXPECTED_REGION_SHAPE), args.compression, args.compression_level
                ),
            )
            regions.attrs["axis_order"] = "event,z,y,x"
            regions.attrs["physical_processing"] = "raw integer counts; no transforms"

            offset = 0
            for schema, keep in zip(schemas, keep_masks, strict=True):
                regions_name = schema.branch_map["regions"]
                assert regions_name is not None
                with uproot.open(schema.path) as root_file:
                    tree = root_file[schema.tree_name]
                    for start in range(0, schema.selected_entries, args.chunk_size):
                        stop = min(start + args.chunk_size, schema.selected_entries)
                        values = np.asarray(
                            tree[regions_name].array(
                                entry_start=start, entry_stop=stop, library="np"
                            )
                        )
                        expected = (stop - start, *EXPECTED_REGION_SHAPE)
                        if values.shape != expected:
                            raise ConversionError(
                                f"Inconsistent chunk shape in {schema.path}:{regions_name}: "
                                f"found {values.shape}, expected {expected}"
                            )
                        if values.dtype.kind not in "iu":
                            raise ConversionError(f"Non-integer raw counts in {schema.path}:{regions_name}")
                        if np.any(values < 0):
                            minimum = int(np.min(values))
                            raise ConversionError(
                                f"Negative raw counts in {schema.path}:{regions_name} (minimum {minimum})"
                            )
                        kept_values = values[keep[start:stop]]
                        rows = kept_values.shape[0]
                        if rows:
                            regions[offset : offset + rows] = kept_values
                            offset += rows
            if offset != total:
                raise ConversionError(f"Wrote {offset} region rows, expected {total}")

            if common_background:
                background = hdf5.create_dataset(
                    "background_mu",
                    data=metadata.background_mu[0].astype(np.float32),
                    **_dataset_kwargs((EXPECTED_REGION_SHAPE[0],), args.compression, args.compression_level),
                )
                background.attrs["common_to_all_events"] = True
            else:
                background = hdf5.create_dataset(
                    "background_mu",
                    data=metadata.background_mu.astype(np.float32),
                    **_dataset_kwargs(
                        (chunk_rows, EXPECTED_REGION_SHAPE[0]),
                        args.compression,
                        args.compression_level,
                    ),
                )
                background.attrs["common_to_all_events"] = False
            background.attrs["axis_order"] = "z" if common_background else "event,z"

            numeric = {
                "presence": metadata.presence,
                "sample_type": metadata.sample_type,
                "slope_xy": metadata.slope_xy,
                "signal_event_id": metadata.signal_event_id,
                "tile_id": metadata.tile_id,
                "crop_id_in_tile": metadata.crop_id_in_tile,
                "negative_type": metadata.negative_type,
                "split": split,
                "source_entry": metadata.source_entry,
            }
            for name, values in numeric.items():
                tail = values.shape[1:]
                dataset = hdf5.create_dataset(
                    name,
                    data=values,
                    **_dataset_kwargs((chunk_rows, *tail), args.compression, args.compression_level),
                )
                if name == "slope_xy":
                    dataset.attrs["component_order"] = "sx,sy"
                    dataset.attrs["unit"] = "xy_bins_per_z_bin"
                elif name == "presence":
                    dataset.attrs["labels"] = "0=no_shower_poisson,1=shower_hard_or_signal"
                elif name == "sample_type":
                    dataset.attrs["mapping_json"] = json.dumps(SAMPLE_TYPE_NAMES)
                elif name == "split":
                    dataset.attrs["mapping_json"] = json.dumps(SPLIT_TO_CODE)
                elif name == "negative_type":
                    dataset.attrs["mapping_json"] = json.dumps(NEGATIVE_TYPE_NAMES)

            string_dtype = h5py.string_dtype(encoding="utf-8")
            hdf5.create_dataset(
                "source_file",
                data=metadata.source_file.astype(string_dtype),
                dtype=string_dtype,
                **_dataset_kwargs((chunk_rows,), args.compression, args.compression_level),
            )

            split_group = hdf5.create_group("split_indices")
            for code, name in enumerate(SPLIT_NAMES):
                indices = np.flatnonzero(split == code).astype(np.int64)
                if indices.size:
                    split_group.create_dataset(
                        name,
                        data=indices,
                        **_dataset_kwargs(
                            (min(chunk_rows, indices.size),), args.compression, args.compression_level
                        ),
                    )
                else:
                    split_group.create_dataset(name, data=indices)

            if metadata.extras:
                extra_group = hdf5.create_group("metadata")
                safe_names: set[str] = set()
                for source_name, values in sorted(metadata.extras.items()):
                    safe_name = _safe_extra_name(source_name)
                    if safe_name in safe_names:
                        raise ConversionError(f"Metadata name collision after HDF5 sanitization: {source_name}")
                    safe_names.add(safe_name)
                    dataset = extra_group.create_dataset(
                        safe_name,
                        data=values,
                        **_dataset_kwargs((chunk_rows,), args.compression, args.compression_level),
                    )
                    dataset.attrs["source_branch"] = source_name

            hdf5.flush()
        os.replace(temporary, output)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return common_background


def write_manifest(path: Path, metadata: Metadata, split: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "dataset_index",
                "split",
                "presence",
                "sample_type",
                "sample_type_name",
                "class_name",
                "group_kind",
                "group_id",
                "signal_event_id",
                "tile_id",
                "crop_id_in_tile",
                "negative_type",
                "negative_type_name",
                "source_file",
                "source_entry",
            ]
        )
        for index in range(len(metadata)):
            sample_type = int(metadata.sample_type[index])
            signal = sample_type == 2
            negative_type = int(metadata.negative_type[index])
            writer.writerow(
                [
                    index,
                    SPLIT_NAMES[int(split[index])],
                    int(metadata.presence[index]),
                    sample_type,
                    SAMPLE_TYPE_NAMES[sample_type],
                    SAMPLE_TYPE_NAMES[sample_type],
                    "signal_event_id" if signal else "tile_id",
                    int(metadata.signal_event_id[index] if signal else metadata.tile_id[index]),
                    int(metadata.signal_event_id[index]),
                    int(metadata.tile_id[index]),
                    int(metadata.crop_id_in_tile[index]),
                    negative_type,
                    NEGATIVE_TYPE_NAMES[negative_type],
                    str(metadata.source_file[index]),
                    int(metadata.source_entry[index]),
                ]
            )


def make_qa_figures(
    hdf5_path: Path, qa_dir: Path, metadata: Metadata, split: np.ndarray, per_class: int
) -> list[str]:
    if per_class <= 0:
        return []
    os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    qa_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []
    with h5py.File(hdf5_path, "r") as hdf5:
        for sample_type, class_name in SAMPLE_TYPE_NAMES.items():
            candidates = np.flatnonzero(metadata.sample_type == sample_type)
            if candidates.size == 0:
                warnings.append(f"No {class_name} samples are available for QA figures")
                continue
            for qa_index, dataset_index in enumerate(candidates[:per_class]):
                volume = hdf5["regions_raw"][dataset_index]
                projection = np.sum(volume, axis=0, dtype=np.int64)
                figure, axis = plt.subplots(figsize=(5.5, 5.0), constrained_layout=True)
                image = axis.imshow(projection, origin="lower", interpolation="nearest", cmap="viridis")
                axis.set_xlabel("x bin")
                axis.set_ylabel("y bin")
                axis.set_title(
                    f"{class_name} #{int(dataset_index)} | {SPLIT_NAMES[int(split[dataset_index])]} | "
                    "raw z-sum"
                )
                figure.colorbar(image, ax=axis, label="raw counts summed over z")
                output = qa_dir / f"{class_name}_{qa_index:02d}_index_{int(dataset_index):06d}.png"
                figure.savefig(output, dpi=150)
                plt.close(figure)
    return warnings


def _schema_payload(schemas: Sequence[RootFileSchema]) -> dict[str, Any]:
    return {
        "expected_hdf5_region_shape": [None, *EXPECTED_REGION_SHAPE],
        "root_axis_order": "entry,z,y,x",
        "files": [
            {
                "path": str(schema.path),
                "tree": schema.tree_name,
                "entries": schema.entries,
                "selected_entries": schema.selected_entries,
                "explicit_split": schema.explicit_split,
                "branch_map": schema.branch_map,
                "branches": {
                    name: {
                        "typename": info.typename,
                        "interpretation": info.interpretation,
                        "dtype": info.dtype,
                        "per_entry_shape": info.per_entry_shape,
                    }
                    for name, info in schema.branches.items()
                },
            }
            for schema in schemas
        ],
    }


def run_inspect(args: argparse.Namespace) -> int:
    schemas = [inspect_root_file(path, args) for path in _expand_inputs(args.inputs)]
    payload = _schema_payload(schemas)
    if args.json:
        print(json.dumps(payload, indent=2, default=_json_default))
    else:
        for schema in schemas:
            print(f"File: {schema.path}")
            print(f"  TTree: {schema.tree_name} ({schema.entries} entries)")
            for name, info in schema.branches.items():
                print(
                    f"  {name}: typename={info.typename}, dtype={info.dtype}, "
                    f"per-entry-shape={info.per_entry_shape}"
                )
            print(f"  resolved roles: {schema.branch_map}")
    return 0


def run_convert(args: argparse.Namespace) -> int:
    schemas = discover_inputs(args.inputs, args)
    metadata, warnings = read_metadata(schemas, args)
    keep_masks, support_filter = build_signal_support_filter(schemas, metadata, args)
    background_tile_filter = apply_background_tile_filter(
        schemas, metadata, keep_masks, args.exclude_background_tile
    )
    keep = np.concatenate(keep_masks)
    metadata = _subset_metadata(metadata, keep)
    if len(metadata) == 0:
        raise ConversionError("No entries remain after the signal-support filter")
    split = assign_splits(metadata, args.seed, args.split_ratios)
    statistics = compute_statistics(metadata, split)
    statistics["signal_support_filter"] = support_filter
    statistics["background_tile_filter"] = background_tile_filter

    output = args.output.resolve()
    manifest = (args.manifest or output.with_suffix(".manifest.csv")).resolve()
    stats_path = (args.statistics or output.with_suffix(".statistics.json")).resolve()
    schema_path = (args.schema_json or output.with_suffix(".schema.json")).resolve()
    qa_dir = (args.qa_dir or output.parent / f"{output.stem}_qa").resolve()

    planned_files = (output, manifest, stats_path, schema_path)
    if len(set(planned_files)) != len(planned_files):
        raise ConversionError("HDF5, manifest, statistics, and schema paths must be distinct")
    source_paths = {schema.path.resolve() for schema in schemas}
    collisions = [str(path) for path in planned_files if path in source_paths]
    if collisions:
        raise ConversionError("Output artifacts cannot overwrite ROOT inputs: " + ", ".join(collisions))
    if not args.overwrite:
        existing = [str(path) for path in planned_files if path.exists()]
        if args.qa_per_class > 0 and qa_dir.exists() and any(qa_dir.iterdir()):
            existing.append(str(qa_dir))
        if existing:
            raise ConversionError(
                "Refusing to overwrite existing conversion artifacts: "
                + ", ".join(existing)
                + "; use --overwrite"
            )
    elif qa_dir.exists():
        for pattern in (
            "signal_*_index_*.png",
            "hard_*_index_*.png",
            "poisson_*_index_*.png",
            # Clean figures produced by schema versions before sample_type became
            # the classification label.
            "positive_*_index_*.png",
            "negative_*_index_*.png",
        ):
            for old_figure in qa_dir.glob(pattern):
                old_figure.unlink()

    common_background = write_hdf5(
        output,
        schemas,
        metadata,
        split,
        keep_masks,
        support_filter,
        background_tile_filter,
        args,
    )
    write_manifest(manifest, metadata, split)
    stats_path.write_text(json.dumps(statistics, indent=2, default=_json_default) + "\n", encoding="utf-8")
    schema_path.write_text(
        json.dumps(_schema_payload(schemas), indent=2, default=_json_default) + "\n",
        encoding="utf-8",
    )
    warnings.extend(make_qa_figures(output, qa_dir, metadata, split, args.qa_per_class))

    total_available = sum(schema.entries for schema in schemas)
    print(f"Converted {len(metadata)} selected entries ({total_available} available in inspected files)")
    print(
        "Excluded unusable signal entries: "
        f"{support_filter['excluded_signal_entries']} "
        f"(minimum support {support_filter['minimum_signal_support_xy']}x"
        f"{support_filter['minimum_signal_support_xy']})"
    )
    print(
        "Excluded background entries from tiles "
        f"{background_tile_filter['excluded_background_tile_ids']}: "
        f"{background_tile_filter['excluded_background_entries']}"
    )
    print(f"HDF5: {output}")
    print(f"Manifest: {manifest}")
    print(f"Statistics: {stats_path}")
    print(f"Schema report: {schema_path}")
    if args.qa_per_class > 0:
        print(f"QA figures: {qa_dir}")
    print(
        "background_mu shape: "
        + (f"{EXPECTED_REGION_SHAPE[:1]} (common to all events)" if common_background else f"({len(metadata)}, 57)")
    )
    print(json.dumps(statistics, indent=2))
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    return 0


def _add_schema_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--tree", default="auto", help="TTree name; default: auto-select the only TTree")
    parser.add_argument("--regions-branch", default="auto")
    parser.add_argument("--background-branch", default="auto")
    parser.add_argument("--presence-branch", default="auto")
    parser.add_argument("--slope-x-branch", default="auto")
    parser.add_argument("--slope-y-branch", default="auto")
    parser.add_argument("--signal-event-id-branch", default="auto")
    parser.add_argument("--tile-id-branch", default="auto")
    parser.add_argument("--crop-id-branch", default="auto")
    parser.add_argument("--sample-type-branch", default="auto")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Inspect and convert raw ROOT source regions to validated HDF5"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect_parser = subparsers.add_parser("inspect", help="show the actual ROOT schema")
    inspect_parser.add_argument("inputs", nargs="+", type=Path)
    inspect_parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    _add_schema_arguments(inspect_parser)
    inspect_parser.set_defaults(handler=run_inspect)

    convert_parser = subparsers.add_parser("convert", help="validate and convert ROOT entries")
    convert_parser.add_argument("inputs", nargs="+", type=Path)
    convert_parser.add_argument("--output", type=Path, required=True)
    convert_parser.add_argument("--manifest", type=Path)
    convert_parser.add_argument("--statistics", type=Path)
    convert_parser.add_argument("--schema-json", type=Path)
    convert_parser.add_argument("--qa-dir", type=Path)
    convert_parser.add_argument("--qa-per-class", type=int, default=3)
    convert_parser.add_argument(
        "--min-signal-support",
        type=int,
        default=20,
        help=(
            "exclude signal entries whose nonzero support is smaller than this in x or y "
            "in any z slice; use 0 to disable (default: 20)"
        ),
    )
    convert_parser.add_argument(
        "--exclude-background-tile",
        "--exclude-hard-tile",
        dest="exclude_background_tile",
        action="append",
        type=int,
        default=list(DEFAULT_EXCLUDED_BACKGROUND_TILES),
        help=(
            "exclude hard and Poisson entries from this tile/cell_id; repeatable "
            "(default: 58; --exclude-hard-tile is a deprecated alias)"
        ),
    )
    convert_parser.add_argument("--max-events", type=int)
    convert_parser.add_argument("--chunk-size", type=int, default=16)
    convert_parser.add_argument("--compression", choices=("gzip", "lzf", "none"), default="gzip")
    convert_parser.add_argument("--compression-level", type=int, default=4)
    convert_parser.add_argument("--split-ratios", nargs=3, type=float, default=(0.70, 0.15, 0.15))
    convert_parser.add_argument("--seed", type=int, default=20260819)
    convert_parser.add_argument("--xy-bin-size-um", type=float, default=50.0)
    convert_parser.add_argument("--z-bin-size-um", type=float, default=1350.0)
    convert_parser.add_argument(
        "--slope-unit-divisor",
        type=float,
        default=1000.0,
        help="source slope divisor to radians (1000 for the inspected mrad branches)",
    )
    convert_parser.add_argument("--overwrite", action="store_true")
    _add_schema_arguments(convert_parser)
    convert_parser.set_defaults(handler=run_convert)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "max_events", None) is not None and args.max_events <= 0:
        parser.error("--max-events must be positive")
    if getattr(args, "chunk_size", 1) <= 0:
        parser.error("--chunk-size must be positive")
    if getattr(args, "qa_per_class", 0) < 0:
        parser.error("--qa-per-class cannot be negative")
    if getattr(args, "min_signal_support", 0) < 0 or getattr(
        args, "min_signal_support", 0
    ) > min(EXPECTED_REGION_SHAPE[1:]):
        parser.error("--min-signal-support must be between 0 and 32")
    if getattr(args, "compression_level", 4) not in range(10):
        parser.error("--compression-level must be between 0 and 9")
    if getattr(args, "xy_bin_size_um", 1.0) <= 0 or getattr(args, "z_bin_size_um", 1.0) <= 0:
        parser.error("bin sizes must be positive")
    if getattr(args, "slope_unit_divisor", 1.0) <= 0:
        parser.error("--slope-unit-divisor must be positive")
    try:
        return int(args.handler(args))
    except ConversionError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
