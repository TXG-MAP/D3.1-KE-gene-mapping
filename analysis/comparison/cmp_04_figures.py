#!/usr/bin/env python3
"""FIG11 — the cross-method comparison figure.

Three panels, in the order the argument runs:

  A  the waterfall, which is what makes the cascade look supported — every step adds Key Events
  B  the conditional test, which refutes it — no step is stronger where the previous ones failed
  C  the mechanism — one shared blind spot, biological level

Panel B is the load-bearing one and is drawn on the Key Events **all three methods attempted**,
not on all 538. Over all 538 the first comparison reverses sign, because it mixes "attempted and
failed" with "never curated" (C40). The restricted denominators are small on the right-hand bar
of each pair, so the counts are printed on the figure rather than left to the caption.

Offline: reads data/comparison/. Run `make comparison` first.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method1"))
import m1lib as m1  # noqa: E402

OUT = m1.ROOT / "data" / "comparison"

LEVEL_ORDER = ["Molecular", "Cellular", "Tissue", "Organ", "Individual", "Population"]
COVERED = "#1baf7a"   # previous step succeeded here
MISSED = "#eb6834"    # previous step failed here — the cascade's target


def panel_a(ax, wf: pd.DataFrame, total: int) -> None:
    labels = ["Method 1", "+ Method 2", "+ Method 3"]
    cum = wf["cumulative"].tolist()
    new = wf["new_kes_added"].tolist()

    prev = 0
    for i, (c, n) in enumerate(zip(cum, new)):
        ax.bar(i, prev, color=m1.GRID, width=0.62, zorder=2)
        ax.bar(i, n, bottom=prev, color=m1.SEQUENTIAL[i + 1], width=0.62, zorder=2)
        ax.text(i, prev + n / 2, f"+{n}", ha="center", va="center",
                color=m1.SURFACE, fontsize=9, fontweight="bold", zorder=3)
        ax.text(i, c + total * 0.022, f"{c}  ({100 * c / total:.0f} %)", ha="center",
                va="bottom", fontsize=8.5, color=m1.INK)
        prev = c

    ax.axhline(total, color=m1.INK_MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)
    ax.text(-0.5, total + total * 0.012, f"all {total} KEs in the four networks", va="bottom",
            ha="left", fontsize=7.8, color=m1.INK_MUTED)
    # Sits in the empty band between the last bar and the 538 line, so the panel does not need
    # a wide right margin — which was squeezing the three tick labels into each other.
    ax.text(2, (cum[-1] + total) / 2 + 6, f"{total - cum[-1]} unreached",
            va="center", ha="center", fontsize=7.8, color=m1.INK_MUTED)

    ax.set_xticks(range(3), labels, fontsize=8)
    ax.set_ylim(0, total * 1.16)
    ax.set_ylabel("Key Events with a gene set")
    ax.set_title("A   Every step adds Key Events", loc="left", fontsize=9.5)
    ax.set_xlim(-0.6, 2.55)
    ax.grid(axis="x", visible=False)
    m1.despine(ax)


def panel_b(ax, tests: pd.DataFrame) -> None:
    rows = tests[tests["universe"] == "KEs all three attempted"]
    rows = rows[rows["previous"].isin(["method1", "method1 or method2"])]
    rows = rows[~((rows["previous"] == "method1") & (rows["next"] == "method3"))]

    labels = ["Method 2\ngiven Method 1", "Method 3\ngiven Methods 1+2"]
    x = range(len(rows))
    w = 0.34
    for i, (_, r) in enumerate(rows.iterrows()):
        ax.bar(i - w / 2, r["b_rate_where_previous_covered"], width=w, color=COVERED, zorder=2)
        ax.bar(i + w / 2, r["b_rate_where_previous_missed"], width=w, color=MISSED, zorder=2)
        ax.text(i - w / 2, r["b_rate_where_previous_covered"] + 2,
                f"{r['b_rate_where_previous_covered']:.0f} %\nn={r['n_previous_covered']}",
                ha="center", va="bottom", fontsize=8, color=m1.INK_SECONDARY,
                linespacing=1.2)
        ax.text(i + w / 2, r["b_rate_where_previous_missed"] + 2,
                f"{r['b_rate_where_previous_missed']:.0f} %\nn={r['n_previous_missed']}",
                ha="center", va="bottom", fontsize=8, color=m1.INK_SECONDARY,
                linespacing=1.2)
        # The cascade predicts the orange bar to be the taller one. It never is.
        ax.annotate("", xy=(i + w / 2, 111), xytext=(i - w / 2, 111),
                    arrowprops=dict(arrowstyle="->", color=m1.INK_MUTED, linewidth=1.0,
                                    connectionstyle="arc3,rad=-0.3"))
        ax.text(i, 116, f"{r['rescue_difference_pp']:+.0f} pp", ha="center", va="bottom",
                fontsize=8.5, color=m1.INK, fontweight="bold")

    ax.set_xticks(list(x), labels, fontsize=8)
    ax.set_ylim(0, 130)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_ylabel("% of Key Events covered")
    ax.set_title("B   but not where the previous one failed", loc="left", fontsize=9.5)
    ax.grid(axis="x", visible=False)
    m1.despine(ax)


def panel_c(ax, lv: pd.DataFrame) -> None:
    lv = lv.reindex([l for l in LEVEL_ORDER if l in lv.index])
    y = range(len(lv))
    ax.barh(list(y), lv["pct_unreached"], color=MISSED, height=0.66, zorder=2)
    for i, (level, r) in enumerate(lv.iterrows()):
        ax.text(r["pct_unreached"] + 1.6, i,
                f"{r['pct_unreached']:.0f} %   ({int(r['unreached'])}/{int(r['kes_in_networks'])})",
                va="center", ha="left", fontsize=8, color=m1.INK_SECONDARY)
    ax.set_yticks(list(y), lv.index)
    ax.invert_yaxis()
    ax.set_xlim(0, 150)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("% of Key Events reached by no method")
    ax.set_title("C   because all three share one blind spot", loc="left", fontsize=9.5)
    ax.grid(axis="y", visible=False)
    m1.despine(ax)


def main() -> None:
    m1.banner("FIG11 — cross-method comparison")
    m1.apply_style()

    wf = pd.read_csv(OUT / "cascade_waterfall.tsv", sep="\t")
    tests = pd.read_csv(OUT / "cascade_rescue_tests.tsv", sep="\t")
    lv = pd.read_csv(OUT / "cascade_residual_by_level.tsv", sep="\t", index_col=0)
    combos = pd.read_csv(OUT / "ke_coverage_combinations.tsv", sep="\t")
    total = int(combos[combos["organ"] == "four-organ union"]["kes"].iloc[0])

    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.5),
                             gridspec_kw={"width_ratios": [1, 1, 1.32], "wspace": 0.5})
    panel_a(axes[0], wf, total)
    panel_b(axes[1], tests)
    panel_c(axes[2], lv)

    fig.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, facecolor=COVERED,
                      label="panel B: where the previous step succeeded"),
        plt.Rectangle((0, 0), 1, 1, facecolor=MISSED,
                      label="where the previous step failed  (what a cascade must rescue)"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.115), ncol=2, fontsize=8.2)

    m1.source_note(fig, -0.215,
                   "Panels A and C cover all 538 Key Events of the four organ networks. Panel "
                   "B is restricted to the 142 Key Events all three methods attempted: over "
                   "all 538 the first comparison reverses sign, because it mixes\nKey Events a "
                   "method attempted and failed with Key Events nobody curated.",
                   fontsize=7.6, linespacing=1.5)

    m1.save_figure(fig, "FIG11-cross-method-comparison")
    plt.close(fig)


if __name__ == "__main__":
    main()
