from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/roc_pr_curves.png"
OUTPUT = ROOT / "output/hmu_review_en/roc_pr_hmu_test_EN.png"


def font(size: int, bold: bool = False):
    candidates = [
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def rounded(draw: ImageDraw.ImageDraw, box, fill, outline):
    draw.rounded_rectangle(box, radius=18, fill=fill, outline=outline, width=2)


def main():
    original = Image.open(SOURCE).convert("RGB")
    width = 1800
    chart_height = round(original.height * width / original.width)
    chart = original.resize((width, chart_height), Image.Resampling.LANCZOS)
    height = chart_height + 250
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    ink = "#102A43"
    blue = "#126782"
    pale = "#EDF3F6"
    muted = "#526777"

    draw.text((72, 44), "H − μ classification performance on test crops", fill=ink, font=font(45, True))
    draw.text(
        (72, 105),
        "ROC and precision–recall curves. Classification includes hard background and signal with true θ > 8 mrad.",
        fill=blue,
        font=font(23),
    )
    rounded(draw, (1180, 38, 1445, 119), pale, blue)
    draw.text((1210, 58), "AUROC  0.99952", fill=ink, font=font(24, True))
    rounded(draw, (1470, 38, 1728, 119), "#F8F4EF", "#C56730")
    draw.text((1497, 58), "AUPRC  0.99924", fill=ink, font=font(24, True))
    draw.line((72, 155, 1728, 155), fill="#C8D8E0", width=2)

    canvas.paste(chart, (0, 170))
    draw.text(
        (72, height - 54),
        "N = 1,730 test crops: 613 signal positives and 1,117 hard-background negatives. Dashed ROC diagonal: random ranking.",
        fill=muted,
        font=font(20),
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(OUTPUT, optimize=True)


if __name__ == "__main__":
    main()
