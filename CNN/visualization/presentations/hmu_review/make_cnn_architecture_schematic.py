from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon


OUT = Path("output/hmu_review_en")
BLUE = "#126782"
INK = "#102A43"
MUTED = "#526777"
PALE = "#EDF3F6"
ORANGE = "#C56730"


def card(
    ax, x, y, w, h, title, body, *, face="#FFFFFF", edge=BLUE,
    title_color=INK, body_size=10.4, compact=False, title_size=None,
):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.012,rounding_size=0.022",
        linewidth=1.8,
        edgecolor=edge,
        facecolor=face,
        zorder=2,
    )
    ax.add_patch(patch)
    title_offset = 0.028 if compact else 0.055
    title_size = title_size if title_size is not None else (9.8 if compact else 12.4)
    body_offset = 0.012 if compact else 0.025
    ax.text(x + w / 2, y + h - title_offset, title, ha="center", va="top", fontsize=title_size,
            color=title_color, fontweight="bold", zorder=3)
    ax.text(x + w / 2, y + h / 2 - body_offset, body, ha="center", va="center", fontsize=body_size,
            color=INK, linespacing=1.25, zorder=3)


def arrow(ax, start, end, *, color=BLUE, width=1.8):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14,
                                 linewidth=width, color=color, zorder=4))


def feature_stack(ax, x, y, w, h, color, layers=4):
    for index in range(layers - 1, -1, -1):
        dx = index * 0.011
        dy = index * 0.010
        ax.add_patch(FancyBboxPatch(
            (x + dx, y + dy), w, h,
            boxstyle="round,pad=0.005,rounding_size=0.014",
            linewidth=0.9,
            edgecolor=color,
            facecolor=color if index == 0 else "#E5F0F3",
            alpha=1.0 if index == 0 else 0.95,
            zorder=2 + (layers - index),
        ))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(16, 9), facecolor="white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.055, 0.948, "Two-Head 2+1D CNN Architecture", fontsize=27, color=INK,
            fontweight="bold", va="top")
    ax.text(0.055, 0.903,
            "H − μ normalized 57 × 20 × 20 crop  →  shower score and direction estimate",
            fontsize=13.5, color=BLUE, va="top")
    ax.plot([0.055, 0.945], [0.875, 0.875], color="#C8D8E0", lw=1.2)

    # Input volume, drawn as stacked XY planes along z.
    ix, iy, iw, ih = 0.065, 0.465, 0.095, 0.235
    for level in range(4):
        offset = level * 0.018
        ax.add_patch(Polygon(
            [[ix + offset, iy + offset], [ix + iw + offset, iy + offset],
             [ix + iw + offset, iy + ih + offset], [ix + offset, iy + ih + offset]],
            closed=True, facecolor="#D7EAF0" if level < 3 else BLUE,
            edgecolor=BLUE, linewidth=1.1, alpha=0.95,
        ))
    # Align both labels with the dark front plane, rather than with the stack origin.
    input_label_x = ix + iw / 2 + 3 * 0.018
    ax.text(input_label_x, iy + ih / 2 + 0.046, "Input", ha="center", va="center",
            color="white", fontsize=12.5, fontweight="bold")
    ax.text(input_label_x, iy + ih / 2 + 0.008, "1 × 57 × 20 × 20", ha="center", va="center",
            color="white", fontsize=7.5)
    ax.text(ix + iw / 2 + 0.028, iy - 0.06, "raw counts after crop\nH − μ normalization", ha="center",
            va="top", color=MUTED, fontsize=10.3, linespacing=1.25)
    arrow(ax, (0.187, 0.595), (0.225, 0.595))

    # Four factorized CNN blocks.
    positions = [0.235, 0.405, 0.575, 0.745]
    channels = [16, 32, 64, 128]
    shapes = ["16 × 57 × 10 × 10", "32 × 57 × 5 × 5", "64 × 28 × 5 × 5", "128 × 14 × 5 × 5"]
    pools = ["MaxPool\n(1, 2, 2)", "MaxPool\n(1, 2, 2)", "MaxPool\n(2, 1, 1)", "MaxPool\n(2, 1, 1)"]
    for index, x in enumerate(positions):
        card(ax, x, 0.445, 0.135, 0.305, f"Block {index + 1} · {channels[index]} channels",
             "Conv 1 × 3 × 3  (XY)\nGroupNorm · SiLU\nConv 3 × 1 × 1  (Z)\nGroupNorm · SiLU\n" + pools[index],
             face="#FFFFFF" if index % 2 else PALE, body_size=9.1, title_size=8.9)
        feature_stack(ax, x + 0.035, 0.365, 0.056, 0.050, BLUE, layers=4)
        ax.text(x + 0.0675, 0.338, shapes[index], ha="center", va="top", fontsize=8.7, color=MUTED)
        if index < 3:
            arrow(ax, (x + 0.135, 0.595), (positions[index + 1] - 0.012, 0.595))

    # Aggregation.
    arrow(ax, (0.882, 0.595), (0.912, 0.595))
    card(ax, 0.912, 0.497, 0.072, 0.195, "Global", "average\n+ max\npooling", face="#F8F4EF", edge=ORANGE,
         title_color=ORANGE, body_size=9.6)
    ax.text(0.948, 0.463, "128 + 128 = 256 features", ha="center", va="top", fontsize=8.9, color=MUTED)

    # Heads.
    ax.text(0.055, 0.290, "Shared 256-feature vector from global pooling", fontsize=12.6, color=BLUE, fontweight="bold")
    ax.text(0.055, 0.245, "Each head: Linear 256 → 64  ·  SiLU  ·  Dropout = 0  ·  Linear", fontsize=11.2,
            color=MUTED)
    ax.add_patch(FancyArrowPatch((0.475, 0.245), (0.535, 0.275), arrowstyle="-|>", mutation_scale=14,
                                 linewidth=1.8, color=ORANGE, connectionstyle="arc3,rad=0.0"))
    ax.add_patch(FancyArrowPatch((0.475, 0.245), (0.535, 0.105), arrowstyle="-|>", mutation_scale=14,
                                 linewidth=1.8, color=ORANGE, connectionstyle="arc3,rad=-0.24"))
    card(ax, 0.55, 0.18, 0.18, 0.17, "Classification head", "1 output\npresence logit → shower score", face="#EEF8FA", edge=BLUE,
         body_size=10.2)
    card(ax, 0.55, 0.03, 0.18, 0.115, "Regression head", "2 outputs: sx, sy", face="#FFF6EE", edge=ORANGE,
         title_color=ORANGE, body_size=9.3, compact=True)
    arrow(ax, (0.73, 0.265), (0.81, 0.265))
    arrow(ax, (0.73, 0.087), (0.81, 0.087), color=ORANGE)
    card(ax, 0.825, 0.185, 0.125, 0.16, "Output", "p(shower)\nclassification", face="#FFFFFF", edge=BLUE,
         body_size=10.2)
    card(ax, 0.825, 0.03, 0.125, 0.115, "Output", "sx, sy", face="#FFFFFF", edge=ORANGE,
         title_color=ORANGE, body_size=9.2, compact=True)

    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"cnn_architecture_schematic_EN.{ext}", dpi=190 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.06, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
