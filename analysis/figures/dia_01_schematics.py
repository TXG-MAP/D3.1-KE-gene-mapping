#!/usr/bin/env python3
"""The two drawn schematics: FIG01 (method overview) and FIG04 (GO BP decision tree).

FIG01 shows the three routes from an AOP-Wiki Key Event to an Ensembl gene set, which part
of the KE record each route reads, and where the three converge on BridgeDb. It is the
orienting figure for Section 3.

FIG04 is the Method 2 decision tree for choosing a GO Biological Process term for a KE's
Process slot. Its logic is the Method 2 selection rubric described in Section 3.2.3 of the
report.

With D31_FIG_CLEAN=1 the standfirst under FIG01 is dropped and both schematics are written
to figures/clean/ instead of figures/final/ — the caption-free variants for the report, which
carries that sentence as the figure caption.

Usage:  python analysis/figures/dia_01_schematics.py
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
FINAL = ROOT / "figures" / "final"
CLEAN = ROOT / "figures" / "clean"

CLEAN_FIGURES = os.environ.get("D31_FIG_CLEAN", "") == "1"

# One palette for both schematics. Method colours are reused in FIG01's lanes.
C = {
    "ke": "#E8EEF7",
    "ke_edge": "#3D5A80",
    "m1": "#FDE2C4",
    "m1_edge": "#C1751B",
    "m2": "#D6ECD9",
    "m2_edge": "#3F7B47",
    "m3": "#E2DAF0",
    "m3_edge": "#5F4B8B",
    "out": "#F2F2F2",
    "out_edge": "#555555",
    "stop": "#F7D9D9",
    "stop_edge": "#B2182B",
    "arrow": "#555555",
    "text": "#1A1A1A",
}


def box(ax, x, y, w, h, text, face, edge, size=9, weight="normal", style="round,pad=0.02"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, linewidth=1.2,
                                facecolor=face, edgecolor=edge, mutation_scale=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size,
            color=C["text"], fontweight=weight, linespacing=1.35, zorder=3)
    return (x, y, w, h)


def arrow(ax, p0, p1, style="-|>", lw=1.3, color=None, connectionstyle=None):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=13,
                                 linewidth=lw, color=color or C["arrow"],
                                 shrinkA=0, shrinkB=0,
                                 connectionstyle=connectionstyle or "arc3,rad=0"))


def save(fig, name):
    target = CLEAN if CLEAN_FIGURES else FINAL
    target.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        out = target / f"{name}.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else None, bbox_inches="tight",
                    facecolor="white")
        print("wrote", out.relative_to(ROOT))
    plt.close(fig)


# --------------------------------------------------------------------------- FIG01

def fig01():
    # The standfirst sits at y=3.2, below the lanes. Dropping it also drops the band it
    # occupied, and the figure height follows so the boxes keep their proportions.
    y0 = 7 if CLEAN_FIGURES else 0
    fig, ax = plt.subplots(figsize=(13, 6.4 * (50 - y0) / 50))
    ax.set_xlim(0, 100)
    ax.set_ylim(y0, 50)
    ax.axis("off")

    # The Key Event record, split into the two parts the methods read.
    box(ax, 0.5, 12, 18, 26, "", C["ke"], C["ke_edge"])
    ax.text(9.5, 36.4, "AOP-Wiki\nKey Event", ha="center", va="top", fontsize=11,
            fontweight="bold", color=C["text"], linespacing=1.25)
    box(ax, 2, 24.5, 15, 5.4, "Free-text\ndescription", "#FFFFFF", C["ke_edge"], size=8.5)
    box(ax, 2, 13.8, 15, 8.8,
        "Key Event Components\n\nProcess · Object · Action",
        "#FFFFFF", C["ke_edge"], size=7.8)

    lanes = [
        (36.5, "Method 1 — text mining", C["m1"], C["m1_edge"],
         "HGNC dictionary\n+ BERN2 NER",
         "Gene symbols named\nin the description"),
        (22.5, "Method 2 — Gene Ontology", C["m2"], C["m2_edge"],
         "GO Biological Process\nterm on the Process slot",
         "Genes annotated\nto the GO term"),
        (8.5, "Method 3 — molecular pathways", C["m3"], C["m3_edge"],
         "molAOP Builder:\nBioBERT ranker + curation",
         "Genes in the\nWikiPathways pathway"),
    ]

    for y, title, face, edge, step, out in lanes:
        ax.text(20.5, y + 9.4, title, ha="left", va="center", fontsize=9.5,
                fontweight="bold", color=edge)
        b1 = box(ax, 20.5, y, 21, 8.2, step, face, edge, size=8.5)
        b2 = box(ax, 46, y, 19, 8.2, out, "#FFFFFF", edge, size=8.5)
        arrow(ax, (b1[0] + b1[2], y + 4.1), (b2[0], y + 4.1), color=edge)
        arrow(ax, (b2[0] + b2[2], y + 4.1), (69.5, y + 4.1), color=edge)

    # Which part of the record each method reads.
    arrow(ax, (17, 27.2), (20.5, 40.6), color=C["m1_edge"],
          connectionstyle="arc3,rad=-0.28")
    arrow(ax, (17, 18.2), (20.5, 26.6), color=C["m2_edge"],
          connectionstyle="arc3,rad=-0.22")
    arrow(ax, (18.5, 33), (20.5, 12.6), color=C["m3_edge"],
          connectionstyle="arc3,rad=0.32")
    ax.text(18.4, 34.6, "description", fontsize=7.2, color=C["m1_edge"], ha="right")
    ax.text(18.4, 22.6, "Process slot", fontsize=7.2, color=C["m2_edge"], ha="right")
    ax.text(20.2, 15.4, "KE title", fontsize=7.2, color=C["m3_edge"], ha="right")

    # Convergence.
    box(ax, 69.5, 8.5, 12, 36.2, "BridgeDb\n\nidentifier\nharmonisation",
        C["out"], C["out_edge"], size=9.5, weight="bold")
    arrow(ax, (81.5, 26.6), (85.5, 26.6))
    box(ax, 85.5, 17, 13, 19.2,
        "Ensembl\ngene set\nper Key Event",
        "#FFFFFF", C["out_edge"], size=9.5, weight="bold")

    if not CLEAN_FIGURES:
        ax.text(50, 3.2,
                "All three methods are run on the same four organ-specific AOP networks "
                "— brain, kidney, liver and lung —\nagainst the pinned AOP-Wiki RDF release, "
                "so any difference between them is a property of the method.",
                ha="center", va="center", fontsize=8.2, color="#444444", linespacing=1.5)

    save(fig, "FIG01-method-overview")


# --------------------------------------------------------------------------- FIG04

def fig04():
    fig, ax = plt.subplots(figsize=(9.2, 12.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 144)
    ax.axis("off")

    steps = [
        (128, "Key Event title, plus the KE and AOP\ndescription where one is available",
         C["ke"], C["ke_edge"], 10, "bold"),
        (117, "Already carries a Biological Process\ncomponent?  If so, keep it.",
         "#FFFFFF", C["ke_edge"], 9, "normal"),
        (105, "Enumerate candidate GO terms\n(AmiGO2 / GOlr search — enumerate, do not recall)",
         C["m2"], C["m2_edge"], 9, "normal"),
        (93, "Rule 1 — drop obsolete classes\nuse the replaced_by term instead",
         "#FFFFFF", C["m2_edge"], 9, "normal"),
        (81, "Rule 2 — keep only directionally neutral terms\n"
             "direction belongs in the PATO Action slot, not in GO\n"
             "unsigned “regulation of X” stays eligible",
         "#FFFFFF", C["m2_edge"], 9, "normal"),
        (68, "Rule 4 — anchor on the KE’s object\n"
             "a KE may take several co-equal Process terms",
         "#FFFFFF", C["m2_edge"], 9, "normal"),
        (56, "Rule 3 — most descriptive, but faithful\n"
             "drop root terms; an exact label match beats a narrower child;\n"
             "do not assert more than the KE does",
         "#FFFFFF", C["m2_edge"], 9, "normal"),
        (43, "Rule 5 — tie-break on relationship tag\nprocess  >  regulation  >  response",
         "#FFFFFF", C["m2_edge"], 9, "normal"),
        (31, "Rule 6 — is the term true of the Key Event?\n"
             "generic is acceptable; untrue is not",
         "#FFFFFF", C["m2_edge"], 9, "normal"),
    ]

    boxes = []
    for y, text, face, edge, size, weight in steps:
        h = 8.4 if text.count("\n") < 2 else 10.2
        boxes.append(box(ax, 12, y, 76, h, text, face, edge, size=size, weight=weight))

    for a, b in zip(boxes, boxes[1:]):
        arrow(ax, (50, a[1]), (50, b[1] + b[3]))

    # The two terminal states.
    box(ax, 4, 15, 40, 9.4,
        "Leave the Process slot empty\n(Rule 6 — a wrong term is worse than none)",
        C["stop"], C["stop_edge"], size=8.5)
    box(ax, 54, 15, 42, 9.4,
        "Assign the GO BP identifier(s)\nto the Key Event",
        C["m2"], C["m2_edge"], size=9, weight="bold")

    last = boxes[-1]
    arrow(ax, (30, last[1]), (24, 24.4))
    arrow(ax, (70, last[1]), (75, 24.4))
    ax.text(23.0, 28.2, "no", fontsize=8.5, color=C["stop_edge"], fontweight="bold")
    ax.text(75.5, 28.2, "yes", fontsize=8.5, color=C["m2_edge"], fontweight="bold")

    box(ax, 54, 3, 42, 8.4,
        "Retrieve the gene set annotated\nto that GO term",
        "#FFFFFF", C["m2_edge"], size=9)
    arrow(ax, (75, 15), (75, 11.4))

    ax.text(50, 140.5,
            "Rule 7 — the KE’s biological level is a soft hint throughout, never a filter",
            ha="center", va="center", fontsize=8.2, color="#444444", style="italic")

    save(fig, "FIG04-go-bp-annotation-decision-tree")


if __name__ == "__main__":
    fig01()
    if CLEAN_FIGURES:
        # FIG04 carries no text below the tree — its only note is the rule-7 hint at the
        # top, which is part of the rubric rather than a caption. Nothing to strip, so no
        # clean variant is written; use figures/final/.
        print("skipping FIG04 — no caption text to strip, use figures/final/")
    else:
        fig04()
