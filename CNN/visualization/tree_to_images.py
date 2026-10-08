#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import numpy as np
import ROOT


def configure_root(palette):
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(palette)
    ROOT.gStyle.SetTitleSize(0, "XYZ")
    ROOT.gStyle.SetLabelSize(0, "XYZ")
    ROOT.gStyle.SetTickLength(0, "XYZ")
    ROOT.gStyle.SetFrameLineWidth(0)
    ROOT.gStyle.SetPadLeftMargin(0)
    ROOT.gStyle.SetPadRightMargin(0)
    ROOT.gStyle.SetPadTopMargin(0)
    ROOT.gStyle.SetPadBottomMargin(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)


def read_sample(input_path, tree_name, branch_name, background_branch, entry):
    root_file = ROOT.TFile.Open(str(input_path))
    tree = root_file.Get(tree_name)
    leaf = tree.GetLeaf(branch_name)
    dimensions = tuple(int(value) for value in re.findall(r"\[(\d+)\]", leaf.GetTitle()))

    tree.GetEntry(entry)
    counts = np.array(getattr(tree, branch_name), dtype=np.int32, copy=True)
    counts = counts.reshape(dimensions)
    background = np.atleast_1d(
        np.array(getattr(tree, background_branch), dtype=np.float32, copy=True)
    )
    if background.size == 1:
        background = np.repeat(background, counts.shape[0])
    root_file.Close()
    return counts, background


def numpy_to_histogram(values, name, bin_size):
    ny, nx = values.shape
    x_half_width = 0.5 * nx * bin_size
    y_half_width = 0.5 * ny * bin_size
    histogram = ROOT.TH2F(
        name,
        "",
        nx,
        -x_half_width,
        x_half_width,
        ny,
        -y_half_width,
        y_half_width,
    )
    histogram.SetDirectory(0)
    for y in range(ny):
        for x in range(nx):
            histogram.SetBinContent(x + 1, y + 1, float(values[y, x]))
    return histogram


def make_smoothed_crops(counts, crop_size, bin_size):
    side_y, side_x = counts.shape[1:]
    first_y = (side_y - crop_size) // 2
    first_x = (side_x - crop_size) // 2
    last_y = first_y + crop_size
    last_x = first_x + crop_size
    result = np.zeros((counts.shape[0], crop_size, crop_size), dtype=np.float32)

    for z in range(counts.shape[0]):
        histogram = numpy_to_histogram(counts[z], f"smooth_source_{z}", bin_size)
        histogram.Smooth()
        for y in range(first_y, last_y):
            for x in range(first_x, last_x):
                result[z, y - first_y, x - first_x] = histogram.GetBinContent(x + 1, y + 1)
    return result


def render_slices(values, background, output_dir, prefix, bin_size, color_max, canvas_size):
    output_dir.mkdir(parents=True, exist_ok=True)
    canvas = ROOT.TCanvas(f"{prefix}_canvas", "", canvas_size, canvas_size)

    for z in range(values.shape[0]):
        histogram = numpy_to_histogram(values[z], f"{prefix}_{z}", bin_size)
        histogram.SetMinimum(float(background[z]))
        histogram.SetMaximum(color_max)
        histogram.Draw("COL0")
        canvas.Update()
        canvas.Print(str(output_dir / f"{prefix}_slice_{z + 1:02d}.png"))
        canvas.Clear()

    canvas.Close()


def write_metadata(output_dir, args, shape, background, raw_max, processed_max):
    metadata = f"""Source ROOT file: {args.input}
Tree: {args.tree}
Entry: {args.entry}
Branch: {args.branch}
Background branch: {args.background_branch}
Branch shape: {list(shape)}
XY bin width: {args.bin_size} micrometres

raw/
- no smoothing
- no crop
- common color minimum: background_mu
- common color maximum: {raw_max}

smooth_crop{args.crop_size}/
- ROOT TH2::Smooth() applied once to each complete source slice
- central crop: {args.crop_size} x {args.crop_size} bins
- common color minimum: background_mu
- common color maximum: {processed_max}

Background values: {background.tolist()}

ROOT palette: {args.palette}
Requested canvas size: {args.canvas_size} x {args.canvas_size} pixels
"""
    (output_dir / "metadata.txt").write_text(metadata)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Render raw and smoothed/cropped images from a ROOT TTree counts branch."
    )
    parser.add_argument("input", type=Path, help="ROOT file containing the samples tree")
    parser.add_argument("--output", type=Path, default=Path("tree_images"))
    parser.add_argument("--tree", default="samples")
    parser.add_argument("--branch", default="counts")
    parser.add_argument("--background-branch", default="background_mu")
    parser.add_argument("--entry", type=int, default=0)
    parser.add_argument("--crop-size", type=int, default=20)
    parser.add_argument("--bin-size", type=float, default=50.0)
    parser.add_argument("--palette", type=int, default=52)
    parser.add_argument("--canvas-size", type=int, default=800)
    return parser.parse_args()


def main():
    args = parse_arguments()
    configure_root(args.palette)
    counts, background = read_sample(
        args.input,
        args.tree,
        args.branch,
        args.background_branch,
        args.entry,
    )

    if counts.ndim != 3:
        raise ValueError(f"Expected a [z,y,x] branch, found shape {counts.shape}")
    if args.crop_size > counts.shape[1] or args.crop_size > counts.shape[2]:
        raise ValueError("crop-size cannot exceed the source slice size")

    processed = make_smoothed_crops(counts, args.crop_size, args.bin_size)
    raw_max = float(np.max(counts))
    processed_max = float(np.max(processed))

    render_slices(
        counts,
        background,
        args.output / "raw",
        "raw",
        args.bin_size,
        raw_max,
        args.canvas_size,
    )
    render_slices(
        processed,
        background,
        args.output / f"smooth_crop{args.crop_size}",
        f"smooth_crop{args.crop_size}",
        args.bin_size,
        processed_max,
        args.canvas_size,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    write_metadata(args.output, args, counts.shape, background, raw_max, processed_max)

    print(f"Input shape: {counts.shape}")
    print(f"Raw color maximum: {raw_max}")
    print(f"Smoothed crop color maximum: {processed_max}")
    print(f"Per-slice color minima from {args.background_branch}: {np.unique(background)}")
    print(f"Images written to: {args.output}")


if __name__ == "__main__":
    main()
