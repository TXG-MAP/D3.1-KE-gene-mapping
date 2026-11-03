#!/usr/bin/env python3
"""Compose the report's screenshot figures from the raw captures.

Raw captures live in figures/source/screenshots/ and are browser-viewport grabs, so each
carries browser furniture and dead space that has to come off before the panels are stacked.
The crop boxes below are the only hand-tuned numbers here; everything downstream (scaling to
a common width, panel letters, borders) is derived.

Builds:
  FIG17  a Key Event in the AOP-Wiki   — free text vs the Key Event Component triple
  FIG19  WikiPathways WP528            — what a Method 3 gene set actually is
  FIG15  the molAOP Builder interface  — candidate ranking, and the admin approval queue

Usage:  python analysis/figures/scr_01_assemble.py
"""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "figures" / "source" / "screenshots"
FINAL = ROOT / "figures" / "final"

OUT_WIDTH = 1500      # common panel width; panels are scaled to it
GAP = 26              # vertical space between panels
MARGIN_L = 60         # left gutter that carries the panel letters
BORDER = "#8A8A8A"
BG = "#FFFFFF"

# (file, crop box, panel letter). Crop boxes are (left, upper, right, lower) on the raw
# capture and trim the browser footer, the dead space below the viewport, and the page
# margins, so each panel is scaled from its own content column.
FIGURES = {
    "FIG17-aopwiki-ke-page": [
        ("FIG17-panelA-ke1539-components.jpg", (0, 0, 1513, 740), "A"),
        ("FIG17-panelB-ke1539-description.jpg", (0, 0, 1513, 740), "B"),
    ],
    "FIG19-wikipathways-wp528": [
        ("FIG19-wikipathways-wp528.jpg", (0, 0, 1513, 672), None),
    ],
    "FIG16-zenodo-deposit": [
        ("FIG16-a-record.jpg", (0, 88, 1513, 636), "A"),
        ("FIG16-b-files.jpg", (0, 0, 1513, 570), "B"),
    ],
    "FIG15-molaop-builder-interface": [
        ("FIG15-a1-ke-selection.jpg", (50, 0, 1140, 440), "A"),
        ("FIG15-a2-suggestions.jpg", (50, 70, 1445, 796), "B"),
        ("FIG15-c-admin-queue.jpg", (25, 0, 1145, 412), "C"),
    ],
}


def _font(size):
    for cand in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    ):
        if Path(cand).exists():
            return ImageFont.truetype(cand, size)
    return ImageFont.load_default()


def build(name, parts):
    panels = []
    for fname, box, letter in parts:
        im = Image.open(SRC / fname).convert("RGB").crop(box)
        scale = OUT_WIDTH / im.width
        im = im.resize((OUT_WIDTH, round(im.height * scale)), Image.LANCZOS)
        panels.append((im, letter))

    total_h = sum(p.height for p, _ in panels) + GAP * (len(panels) - 1)
    canvas = Image.new("RGB", (MARGIN_L + OUT_WIDTH + 8, total_h + 8), BG)
    draw = ImageDraw.Draw(canvas)
    font = _font(34)

    y = 4
    for im, letter in panels:
        canvas.paste(im, (MARGIN_L, y))
        draw.rectangle([MARGIN_L - 1, y - 1, MARGIN_L + OUT_WIDTH, y + im.height],
                       outline=BORDER, width=1)
        if letter:
            draw.text((12, y + 2), letter, fill="#000000", font=font)
        y += im.height + GAP

    png = FINAL / f"{name}.png"
    canvas.save(png, dpi=(300, 300))
    canvas.save(FINAL / f"{name}.pdf", "PDF", resolution=300)
    print(f"wrote {png.relative_to(ROOT)}  ({canvas.width}x{canvas.height})")


def main():
    FINAL.mkdir(parents=True, exist_ok=True)
    for name, parts in FIGURES.items():
        build(name, parts)


if __name__ == "__main__":
    main()
