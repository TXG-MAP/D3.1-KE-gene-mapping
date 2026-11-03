#!/usr/bin/env python3
"""Stage 5 — the Method 1 figures.

Reads only the derived tables from stages 1-4, so figures can be restyled without touching
the endpoint.

Design constraints applied throughout:

* Four fixed organ hues, assigned by organ and never cycled or reassigned when a subset is
  drawn. The set validates under ``--pairs all`` (worst CVD deltaE 9.2, worst normal-vision
  deltaE 16.3 against the light surface).
* The aqua used for liver sits just below 3:1 contrast on the surface, so the documented
  relief applies: every figure that uses it carries direct value labels, and each figure is
  backed by a CSV table.
* One axis per panel — no dual scales anywhere.
* Recessive grid and axes; values and labels in ink tokens, never in the series colour.
* A legend whenever two or more series are drawn.

Usage:
    python m1_05_figures.py
"""

from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse, Patch

import m1lib as m1

m1.apply_style()

STATS = json.loads((m1.DERIVED_DIR / "method1_statistics.json").read_text())


def _pairs() -> pd.DataFrame:
    return pd.read_csv(m1.RAW_DIR / "ke_gene_long.csv", dtype={"hgnc": str})


def _frame_kes(frame: str, organ: str) -> set[str]:
    path = m1.FRAME_DIR / f"{frame}_{organ.lower()}.csv"
    return set(pd.read_csv(path)["keid"]) if path.exists() else set()


def _label_bars(ax, bars, values, fmt="{:.0f}", inside_threshold=0.18, horizontal=True):
    """Direct value labels — the contrast relief the palette requires, and the thing that
    lets a reader use the chart without going back to the table."""
    span = max(values) if len(values) else 1
    for bar, value in zip(bars, values):
        if horizontal:
            inside = value / span > inside_threshold
            x = bar.get_width() - span * 0.015 if inside else bar.get_width() + span * 0.015
            ax.text(x, bar.get_y() + bar.get_height() / 2, fmt.format(value),
                    va="center", ha="right" if inside else "left", fontsize=8,
                    color=m1.SURFACE if inside else m1.INK_SECONDARY,
                    fontweight="bold" if inside else "normal")
        else:
            inside = value / span > inside_threshold
            y = bar.get_height() - span * 0.02 if inside else bar.get_height() + span * 0.02
            ax.text(bar.get_x() + bar.get_width() / 2, y, fmt.format(value),
                    ha="center", va="top" if inside else "bottom", fontsize=8,
                    color=m1.SURFACE if inside else m1.INK_SECONDARY,
                    fontweight="bold" if inside else "normal")


# --------------------------------------------------------------------------------------
# FIG07 — four-organ gene overlap
# --------------------------------------------------------------------------------------

# Four-ellipse Venn geometry. Four sets cannot be drawn with circles, and the ellipse layout
# below is the standard construction; region label positions are keyed by a membership
# bitstring over the four sets in order.
_ELLIPSES = [(0.350, 0.400, 0.72, 0.45, 140.0),
             (0.450, 0.500, 0.72, 0.45, 140.0),
             (0.544, 0.500, 0.72, 0.45, 40.0),
             (0.644, 0.400, 0.72, 0.45, 40.0)]
_REGION_XY = {
    "0001": (0.87, 0.42), "0010": (0.68, 0.72), "0011": (0.77, 0.59), "0100": (0.32, 0.72),
    "0101": (0.71, 0.30), "0110": (0.50, 0.66), "0111": (0.65, 0.50), "1000": (0.14, 0.42),
    "1001": (0.50, 0.17), "1010": (0.29, 0.30), "1011": (0.39, 0.24), "1100": (0.23, 0.59),
    "1101": (0.61, 0.24), "1110": (0.35, 0.50), "1111": (0.50, 0.38),
}
_SET_LABEL_XY = [(0.13, 0.18), (0.18, 0.83), (0.82, 0.83), (0.87, 0.18)]


def fig07_overlap() -> None:
    pairs = _pairs()
    sets = {o: set(pairs.loc[pairs.keid.isin(_frame_kes("sidra", o)), "hgnc"])
            for o in m1.ORGANS}
    organs = list(m1.ORGANS)

    fig = plt.figure(figsize=(11.5, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.22)

    # --- panel A: the Venn -------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.grid(False)
    for (x, y, w, h, angle), organ in zip(_ELLIPSES, organs):
        colour = m1.ORGAN_COLOURS[organ]
        ax.add_patch(Ellipse((x, y), w, h, angle=angle, facecolor=colour, alpha=0.22,
                             edgecolor="none"))
        ax.add_patch(Ellipse((x, y), w, h, angle=angle, facecolor="none", edgecolor=colour,
                             linewidth=1.6))
    for bits, (x, y) in _REGION_XY.items():
        members = [organs[i] for i, b in enumerate(bits) if b == "1"]
        inside = set.intersection(*(sets[o] for o in members))
        outside = set().union(*(sets[o] for o in organs if o not in members)) \
            if len(members) < 4 else set()
        n = len(inside - outside)
        if n:
            ax.text(x, y, str(n), ha="center", va="center", fontsize=9, color=m1.INK,
                    fontweight="bold" if len(members) == 4 else "normal")
    for (x, y), organ in zip(_SET_LABEL_XY, organs):
        ax.text(x, y, organ, ha="center", va="center", fontsize=10, fontweight="bold",
                color=m1.ORGAN_COLOURS[organ])
    ax.set_title("A  Genes shared between the organ networks", loc="left")

    # --- panel B: the same regions, sorted, where they can actually be read -------------
    ax2 = fig.add_subplot(gs[0, 1])
    regions = []
    for bits in _REGION_XY:
        members = [organs[i] for i, b in enumerate(bits) if b == "1"]
        inside = set.intersection(*(sets[o] for o in members))
        outside = set().union(*(sets[o] for o in organs if o not in members)) \
            if len(members) < 4 else set()
        n = len(inside - outside)
        if n:
            regions.append((" & ".join(members), n, members))
    regions.sort(key=lambda r: r[1])
    labels = [r[0] for r in regions]
    values = [r[1] for r in regions]
    # A region belonging to one organ takes that organ's hue; a shared region is not any one
    # organ's, so it takes a neutral ink rather than an arbitrary member's colour.
    colours = [m1.ORGAN_COLOURS[r[2][0]] if len(r[2]) == 1 else m1.INK_MUTED for r in regions]
    bars = ax2.barh(range(len(regions)), values, color=colours, height=0.72)
    _label_bars(ax2, bars, values)
    ax2.set_yticks(range(len(regions)), labels, fontsize=8)
    ax2.set_xlabel("genes in this region and no other")
    ax2.set_xlim(0, max(values) * 1.12)
    ax2.grid(axis="y", visible=False)
    m1.despine(ax2, keep=("bottom",))
    ax2.set_title("B  Region sizes, sorted", loc="left")
    ax2.legend(handles=[Patch(facecolor=m1.INK_MUTED, label="shared between networks")] +
                       [Patch(facecolor=m1.ORGAN_COLOURS[o], label=f"{o} only")
                        for o in organs],
               loc="lower right", ncol=1, fontsize=8)

    fig.suptitle("Overlap of Method 1 gene sets across the four organ AOP networks",
                 x=0.008, ha="left", fontsize=12, fontweight="bold")
    m1.source_note(fig, -0.02,
                   "AOP-Wiki RDF 97191db (2026-08-01), regex + BERN2 union, provisional "
                   "(February 2024) organ frame.")
    m1.save_figure(fig, "FIG07-method1-gene-venn-four-organs")
    plt.close(fig)


# --------------------------------------------------------------------------------------
# FIG08 — regex against BERN2
# --------------------------------------------------------------------------------------


def fig08_detectors() -> None:
    decomposition = pd.DataFrame(STATS["detector_decomposition"])
    split = pd.DataFrame(STATS["organ_detector_split"])
    split = split[split.frame == "sidra"].set_index("organ").reindex(m1.ORGANS)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw={"wspace": 0.26})

    # --- panel A: repository-wide, the three ways of counting it ------------------------
    # Shown as shares, not counts. The three rows have totals of 488, 1914 and 4366, so on a
    # shared count axis the Key Event row collapses to a sliver and the panel answers "which
    # row is biggest" — which is not a question anyone has. The question is what share each
    # detector contributed, so the bars are normalised and the counts printed on the segments.
    ax = axes[0]
    order = ["regex only", "both", "BERN2 only"]
    y = np.arange(len(decomposition))
    totals = decomposition[order].sum(axis=1).values
    left = np.zeros(len(decomposition))
    for key in order:
        counts = decomposition[key].values
        shares = 100 * counts / totals
        bars = ax.barh(y, shares, left=left, height=0.6, color=m1.DETECTOR_COLOURS[key],
                       edgecolor=m1.SURFACE, linewidth=2)  # 2px surface gap between segments
        for bar, share, count, base in zip(bars, shares, counts, left):
            if share > 6:
                ax.text(base + share / 2, bar.get_y() + bar.get_height() / 2,
                        f"{count:,}\n{share:.0f} %", ha="center", va="center", fontsize=8,
                        color=m1.SURFACE, fontweight="bold", linespacing=1.25)
        left += shares
    ax.set_yticks(y, [f"{label}\n(n = {int(total):,})"
                      for label, total in zip(decomposition["counted over"], totals)])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of the union, %")
    ax.grid(axis="y", visible=False)
    m1.despine(ax, keep=("bottom",))
    ax.set_title("A  Repository-wide, counted three ways", loc="left")
    ax.legend(handles=[Patch(facecolor=m1.DETECTOR_COLOURS[k], label=k) for k in order],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=8.5)

    # --- panel B: per organ network -----------------------------------------------------
    ax2 = axes[1]
    y = np.arange(len(split))
    left = np.zeros(len(split))
    totals = split[order].sum(axis=1).values
    for key in order:
        widths = split[key].values
        bars = ax2.barh(y, widths, left=left, height=0.6, color=m1.DETECTOR_COLOURS[key],
                        edgecolor=m1.SURFACE, linewidth=2)
        for bar, width, base in zip(bars, widths, left):
            if width / totals.max() > 0.06:
                ax2.text(base + width / 2, bar.get_y() + bar.get_height() / 2, f"{width:,}",
                         ha="center", va="center", fontsize=8.5, color=m1.SURFACE,
                         fontweight="bold")
        left += widths
    for i, total in enumerate(totals):
        ax2.text(total + totals.max() * 0.015, i, f"{total:,}", va="center", fontsize=8,
                 color=m1.INK_SECONDARY)
    ax2.set_yticks(y, split.index)
    ax2.invert_yaxis()
    ax2.set_xlim(0, totals.max() * 1.1)
    ax2.set_xlabel("unique genes")
    ax2.grid(axis="y", visible=False)
    m1.despine(ax2, keep=("bottom",))
    ax2.set_title("B  Per organ network", loc="left")

    fig.suptitle("Genes found by the HGNC dictionary, by BERN2, and by both",
                 x=0.008, ha="left", fontsize=12, fontweight="bold")
    m1.source_note(fig, -0.13,
                   "AOP-Wiki RDF 97191db (2026-08-01). Panel B uses the provisional "
                   "(February 2024) organ frame; totals are unique genes per network.")
    m1.save_figure(fig, "FIG08-method1-regex-vs-bern2")
    plt.close(fig)


# --------------------------------------------------------------------------------------
# FIG09 — genes per Key Event
# --------------------------------------------------------------------------------------


def fig09_genes_per_ke() -> None:
    pairs = _pairs()
    data, labels, colours = [], [], []
    for organ in m1.ORGANS:
        counts = (pairs[pairs.keid.isin(_frame_kes("sidra", organ))]
                  .groupby("keid")["hgnc"].nunique().values)
        data.append(counts)
        labels.append(organ)
        colours.append(m1.ORGAN_COLOURS[organ])

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    positions = np.arange(len(data))
    box = ax.boxplot(data, positions=positions, widths=0.52, patch_artist=True,
                     showfliers=False, medianprops=dict(color=m1.INK, linewidth=2),
                     whiskerprops=dict(color=m1.INK_MUTED, linewidth=1),
                     capprops=dict(color=m1.INK_MUTED, linewidth=1))
    for patch, colour in zip(box["boxes"], colours):
        patch.set_facecolor(colour)
        patch.set_alpha(0.28)
        patch.set_edgecolor(colour)
        patch.set_linewidth(1.6)

    rng = np.random.default_rng(20260803)  # fixed seed: the jitter must not move between runs
    for pos, counts, colour in zip(positions, data, colours):
        x = pos + rng.uniform(-0.16, 0.16, size=len(counts))
        ax.scatter(x, counts, s=13, color=colour, alpha=0.55, linewidths=0.5,
                   edgecolors=m1.SURFACE, zorder=3)

    for pos, counts in zip(positions, data):
        median = float(np.median(counts))
        ax.text(pos, ax.get_ylim()[1] * 0.97, f"median {median:.0f}\nn = {len(counts)}",
                ha="center", va="top", fontsize=8.5, color=m1.INK_SECONDARY)

    ax.set_xticks(positions, labels)
    ax.set_ylabel("genes per annotated Key Event")
    ax.grid(axis="x", visible=False)
    m1.despine(ax)
    kw = STATS["organ_tests"]["genes_per_ke_kruskal"]
    ax.set_title("Genes per Key Event, by organ AOP network", loc="left")
    m1.source_note(fig, -0.03,
                   f"Kruskal-Wallis across the four networks: H = {kw['H']:.2f}, "
                   f"p = {kw['p_value']:.2f}; no pairwise difference\nsurvives Holm "
                   "correction. Only Key Events carrying at least one gene are shown.\n"
                   "Outliers are drawn as points rather than repeated as fliers.")
    m1.save_figure(fig, "FIG09-genes-per-ke-boxplot")
    plt.close(fig)


# --------------------------------------------------------------------------------------
# FIG12 — where the method works and where it does not
# --------------------------------------------------------------------------------------


def fig12_structure() -> None:
    """Coverage by structure, against BOTH denominators.

    Two series per category rather than one, because the obvious objection to the raw gradient
    is that molecular Key Events might simply be better described. Plotting only the "all KEs"
    series would leave that unanswered; plotting only the conditional series would hide how
    much of the shortfall is missing text rather than method failure.
    """
    lvl = pd.DataFrame(STATS["coverage_by_biolevel"]["table"])
    role = pd.DataFrame(STATS["coverage_by_role"]["table"])
    chi = STATS["coverage_by_biolevel"]

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8),
                             gridspec_kw={"wspace": 0.3, "width_ratios": [1.3, 1.0]})

    for ax, df, cat_col, title in (
        (axes[0], lvl, "biological level", "A  By biological level of organisation"),
        (axes[1], role, "role", "B  By role in the AOP"),
    ):
        y = np.arange(len(df))
        height = 0.36
        cond = df["% annotated of described"].values
        allk = df["% annotated of all"].values if "% annotated of all" in df \
            else df["% annotated"].values
        b1 = ax.barh(y - height / 2, cond, height=height, color=m1.SEQUENTIAL[3])
        b2 = ax.barh(y + height / 2, allk, height=height, color=m1.SEQUENTIAL[0])
        _label_bars(ax, b1, cond, fmt="{:.0f} %")
        _label_bars(ax, b2, allk, fmt="{:.0f} %")
        ax.set_yticks(y, [f"{name}\nn = {int(n):,}"
                          for name, n in zip(df[cat_col], df["KEs"])])
        ax.invert_yaxis()
        ax.set_xlim(0, 105)
        ax.set_xlabel("% of Key Events carrying at least one gene")
        ax.grid(axis="y", visible=False)
        m1.despine(ax, keep=("bottom",))
        ax.set_title(title, loc="left")

    axes[0].legend(handles=[
        Patch(facecolor=m1.SEQUENTIAL[3], label="of Key Events that have a description"),
        Patch(facecolor=m1.SEQUENTIAL[0], label="of all Key Events"),
    ], loc="upper center", bbox_to_anchor=(0.75, -0.20), ncol=2, fontsize=8.5)

    fig.suptitle("Method 1 finds genes where the text is about molecules",
                 x=0.008, ha="left", fontsize=12, fontweight="bold")
    m1.source_note(fig, -0.26,
                   f"All 1595 Key Events of AOP-Wiki RDF 97191db.\n"
                   f"A, among described Key Events: chi2({chi['conditional_dof']}) = "
                   f"{chi['conditional_chi2']:.0f}, p = {chi['conditional_p_value']:.1e} "
                   f"(over all KEs: chi2 = {chi['chi2']:.0f}, p = {chi['p_value']:.1e}).\n"
                   "Conditioning on having a description does not explain the gradient away "
                   "— it roughly doubles the test statistic.\n"
                   "Roles in B overlap, so those three bars are not a partition of the 1595.")
    m1.save_figure(fig, "FIG12-method1-coverage-by-structure")
    plt.close(fig)


# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# FIG14 — the denominator
# --------------------------------------------------------------------------------------


def fig14_denominator() -> None:
    d = STATS["denominator"]
    quart = pd.DataFrame(STATS["description_length"]["quartiles"])

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.0),
                             gridspec_kw={"wspace": 0.28, "width_ratios": [1.1, 1.0]})

    # --- panel A: what the 1595 actually decomposes into --------------------------------
    ax = axes[0]
    segments = [
        ("annotated by Method 1", d["annotated"], m1.SEQUENTIAL[3]),
        ("described, no gene found", d["kes_with_description"] - d["annotated"],
         m1.SEQUENTIAL[1]),
        ("no description to mine", d["kes_without_description"], m1.INK_MUTED),
    ]
    left = 0
    for label, value, colour in segments:
        ax.barh([0], [value], left=[left], height=0.6, color=colour,
                edgecolor=m1.SURFACE, linewidth=2)
        ax.text(left + value / 2, 0, f"{value:,}\n{100 * value / d['kes_total']:.0f} %",
                ha="center", va="center", fontsize=9, color=m1.SURFACE, fontweight="bold",
                linespacing=1.3)
        left += value
    # The two denominators, drawn as spans beneath the bar so the reader can see that the
    # same 488 sits over two different bases.
    ax.plot([0, d["kes_total"]], [-0.46, -0.46], color=m1.INK_SECONDARY, linewidth=1.5)
    ax.text(d["kes_total"] / 2, -0.54,
            f"all Key Events (n = {d['kes_total']:,})   ->   "
            f"{100 * d['coverage_all_kes']:.1f} % coverage",
            ha="center", va="top", fontsize=8.5, color=m1.INK_SECONDARY)
    ax.plot([0, d["kes_with_description"]], [0.46, 0.46], color=m1.INK, linewidth=1.5)
    ax.text(d["kes_with_description"] / 2, 0.54,
            f"Key Events the method can read (n = {d['kes_with_description']:,})   ->   "
            f"{100 * d['coverage_described_kes']:.1f} % coverage",
            ha="center", va="bottom", fontsize=8.5, color=m1.INK, fontweight="bold")
    ax.set_ylim(-0.85, 0.85)
    ax.set_xlim(0, d["kes_total"] * 1.01)
    ax.set_yticks([])
    ax.set_xlabel("Key Events")
    ax.grid(visible=False)
    m1.despine(ax, keep=("bottom",))
    ax.set_title("A  The same 488 Key Events, two denominators", loc="left")
    ax.legend(handles=[Patch(facecolor=c, label=l) for l, _, c in segments],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=8.5)

    # --- panel B: how much text it takes ------------------------------------------------
    ax2 = axes[1]
    x = np.arange(len(quart))
    bars = ax2.bar(x, quart["% annotated"], width=0.62, color=m1.SEQUENTIAL[2])
    _label_bars(ax2, bars, quart["% annotated"].values, fmt="{:.0f} %", horizontal=False)
    for i, row in quart.iterrows():
        ax2.text(i, 3, f"median\n{int(row['median_genes'])} genes", ha="center", va="bottom",
                 fontsize=8, color=m1.SURFACE)
    ax2.set_xticks(x, [f"{q}\n~{int(l):,} chars"
                       for q, l in zip(quart["description length"], quart["median_length"])],
                   fontsize=8)
    ax2.set_ylim(0, 108)
    ax2.set_ylabel("% of Key Events carrying at least one gene")
    ax2.grid(axis="x", visible=False)
    m1.despine(ax2)
    ax2.set_title("B  Longer descriptions, more genes", loc="left")

    dl = STATS["description_length"]
    fig.suptitle("Method 1's ceiling is how much text AOP-Wiki has, not the method",
                 x=0.008, ha="left", fontsize=12, fontweight="bold")
    m1.source_note(fig, -0.52,
                   f"AOP-Wiki RDF 97191db (2026-08-01). Every gene-annotated Key Event has a "
                   f"description and no undescribed one carries a gene,\nso the "
                   f"{d['kes_without_description']:,} Key Events without a description are "
                   f"out of scope rather than missed.\nB: description length against genes "
                   f"found, Spearman rho = {dl['spearman_rho']:.2f} "
                   f"(p = {dl['spearman_p']:.0e}); "
                   f"quartiles of the {d['kes_with_description']} described Key Events.")
    m1.save_figure(fig, "FIG14-method1-denominator")
    plt.close(fig)




def main() -> None:
    m1.banner("figures")
    fig07_overlap()
    fig08_detectors()
    fig09_genes_per_ke()
    fig12_structure()
    fig14_denominator()
    # The manifest is the provenance record for figures/final/. A caption-free run writes
    # nowhere near it, so it must not restamp it either.
    if not m1.CLEAN_FIGURES:
        m1.write_manifest("05_figures", {"figures": [
            "FIG07-method1-gene-venn-four-organs", "FIG08-method1-regex-vs-bern2",
            "FIG09-genes-per-ke-boxplot", "FIG12-method1-coverage-by-structure",
            "FIG14-method1-denominator",
        ]})
    print("\nStage 5 complete.")


if __name__ == "__main__":
    main()
