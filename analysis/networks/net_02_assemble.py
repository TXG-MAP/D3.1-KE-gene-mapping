#!/usr/bin/env python3
"""FIG18 — compose the four Cytoscape renderings into the report figure.

Reads the per-organ PNGs written by net_01_cytoscape_render.py, trims their whitespace,
lays them out as a 2x2 panel with a shared role legend, and writes
figures/final/FIG18-organ-network-graph.{png,pdf}.

Panel titles carry the KE and KER counts read back from data/networks/, so the figure
cannot drift from Table 1 without the numbers here changing too.

Usage:  python analysis/networks/net_02_assemble.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
NETDIR = ROOT / "data" / "networks"
SRC = ROOT / "figures" / "source" / "networks"
FINAL = ROOT / "figures" / "final"

# Panel order and the organ names as the report uses them. The working files say "lungs";
# the report says "lung" in prose.
PANELS = [
    ("A", "brain", "Brain"),
    ("B", "kidney", "Kidney"),
    ("C", "liver", "Liver"),
    ("D", "lungs", "Lung"),
]

ROLE_COLOURS = {
    "Molecular Initiating Event": "#1B7837",
    "Key Event": "#BFD3E6",
    "Adverse Outcome": "#B2182B",
}


def trim(path, pad=12):
    """Crop the uniform white margin Cytoscape leaves around the drawn network."""
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    bbox = Image.new("RGB", im.size)
    from PIL import ImageChops

    diff = ImageChops.difference(im, bg)
    box = diff.getbbox()
    if box is None:
        return im
    left, upper, right, lower = box
    left, upper = max(0, left - pad), max(0, upper - pad)
    right, lower = min(im.width, right + pad), min(im.height, lower + pad)
    return im.crop((left, upper, right, lower))


def counts(organ):
    nodes = pd.read_csv(NETDIR / f"{organ}_nodes.tsv", sep="\t", dtype=str)
    edges = pd.read_csv(NETDIR / f"{organ}_edges.tsv", sep="\t", dtype=str)
    n = nodes["type"].value_counts()
    return len(nodes), len(edges), int(n.get("MIE", 0)), int(n.get("AO", 0))


def main():
    FINAL.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.08,
                        wspace=0.04, hspace=0.14)

    for ax, (letter, organ, label) in zip(axes.ravel(), PANELS):
        im = trim(SRC / f"{organ}.png")
        ax.imshow(im)
        ax.set_axis_off()

        nk, nr, mie, ao = counts(organ)
        ax.set_title(
            f"{label}  —  {nk} KEs, {nr} KERs  ({mie} MIEs, {ao} AOs)",
            fontsize=10.5, pad=10,
        )
        ax.text(-0.02, 1.055, letter, transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="bottom", ha="left")

    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=10,
               markerfacecolor=c, markeredgecolor="#444444", label=k)
        for k, c in ROLE_COLOURS.items()
    ]
    handles.append(
        Line2D([], [], color="#9E9E9E", linewidth=1.4,
               label="Key Event Relationship (direction of causation)")
    )
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9.5, handletextpad=0.5, columnspacing=1.6,
               bbox_to_anchor=(0.5, 0.005))

    for ext in ("png", "pdf"):
        out = FINAL / f"FIG18-organ-network-graph.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else None, bbox_inches="tight",
                    facecolor="white")
        print("wrote", out.relative_to(ROOT))


if __name__ == "__main__":
    main()
