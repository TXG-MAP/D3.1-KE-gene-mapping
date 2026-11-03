#!/usr/bin/env python3
"""FIG10 — what Method 3 delivered: the approved KE→WikiPathways mappings.

Two panels, both computed from the mappings the molAOP Builder actually holds:

  A  coverage of each organ network
  B  coverage by biological level, which is where the method's reach runs out

**Scope decision (Marvin, 2026-08-14): this figure reports the DELIVERED mappings, not the
statistics of the adjudication run that produced some of them.** An earlier version of FIG10
showed the funnel from ranked suggestion to approved mapping and the reasons 3914 candidates
were declined. Those are properties of one curation exercise rather than of the resource D3.1
reports on, and they are out of scope here. The underlying records are kept in
data/method3/derived/assessed_candidates.tsv for D3.2, which assesses the completeness of this
method and is where they belong.

Panel B is the inverse view of FIG11 panel C and is deliberately kept: FIG11C shows what NO
method reaches across all three, this shows how far Method 3 alone gets at each level.

Offline: reads data/method3/derived/. Run `make method3-coverage` first.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method1"))
import m1lib as m1  # noqa: E402

ROOT = m1.ROOT
DERIVED = ROOT / "data" / "method3" / "derived"

LEVELS = ["Molecular", "Cellular", "Tissue", "Organ", "Individual", "Population"]
NODE_FILE = {"brain": "brain_nodes.tsv", "kidney": "kidney_nodes.tsv",
             "liver": "liver_nodes.tsv", "lung": "lungs_nodes.tsv"}


def panel_a(ax, cov: pd.DataFrame) -> None:
    rows = cov[cov["organ"] != "four-organ union"].copy()
    union = cov[cov["organ"] == "four-organ union"].iloc[0]

    labels = [o.capitalize() for o in rows["organ"]] + ["All four\nnetworks"]
    values = list(rows["pct_kes_covered"]) + [union["pct_kes_covered"]]
    counts = list(zip(rows["kes_with_approved_mapping"], rows["kes_in_frame"])) + \
        [(union["kes_with_approved_mapping"], union["kes_in_frame"])]
    colours = [m1.ORGAN_COLOURS[o.capitalize()] for o in rows["organ"]] + [m1.INK_SECONDARY]

    x = range(len(values))
    ax.bar(list(x), values, color=colours, width=0.62, zorder=2)
    for i, (v, (n, d)) in enumerate(zip(values, counts)):
        ax.text(i, v + 1.6, f"{v:.1f} %", ha="center", va="bottom", fontsize=8.5,
                color=m1.INK, fontweight="bold" if i == len(values) - 1 else "normal")
        ax.text(i, v / 2, f"{n}\nof {d}", ha="center", va="center", fontsize=8,
                color=m1.SURFACE, linespacing=1.3)

    ax.set_xticks(list(x), labels, fontsize=8.5)
    ax.set_ylim(0, 78)
    ax.set_ylabel("% of the network's Key Events\nwith an approved mapping")
    ax.set_title("A   Coverage of the four organ networks", loc="left", fontsize=9.5)
    ax.grid(axis="x", visible=False)
    m1.despine(ax)


def panel_b(ax, lv: pd.DataFrame) -> None:
    y = range(len(lv))
    ax.barh(list(y), lv["pct"], color=m1.SEQUENTIAL[2], height=0.66, zorder=2)
    for i, r in enumerate(lv.itertuples()):
        ax.text(r.pct + 1.5, i, f"{r.pct:.1f} %   ({r.mapped}/{r.kes})", va="center",
                ha="left", fontsize=8, color=m1.INK_SECONDARY)
    ax.set_yticks(list(y), lv["level"], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 106)
    ax.set_xticks([0, 25, 50, 75])
    ax.set_xlabel("% of Key Events with an approved mapping")
    ax.set_title("B   and where its reach runs out", loc="left", fontsize=9.5)
    ax.grid(axis="y", visible=False)
    m1.despine(ax)


def main() -> None:
    m1.banner("FIG10 — Method 3 delivered mappings")
    m1.apply_style()

    cov = pd.read_csv(DERIVED / "method3_coverage_by_organ.tsv", sep="\t")
    pairs = pd.read_csv(DERIVED / "method3_approved_pairs.tsv", sep="\t")

    networks = {o: set(pd.read_csv(ROOT / "data" / "networks" / f, sep="\t")["id"].str.strip())
                for o, f in NODE_FILE.items()}
    union = set().union(*networks.values())
    mapped = set(pairs[pairs["ke_id"].isin(union)]["ke_id"])

    # Biological level comes from the release inventory, not from the builder: the builder's
    # own biolevel string is free text and varies in case and wording between records.
    inv = pd.read_csv(ROOT / "data/method1/raw/ke_inventory.csv")
    inv["keid"] = inv["keid"].str.strip()
    level_of = dict(zip(inv["keid"], inv["biolevel"]))

    rows = []
    for level in LEVELS:
        kes = {k for k in union if level_of.get(k) == level}
        rows.append({"level": level, "kes": len(kes), "mapped": len(kes & mapped),
                     "pct": round(100 * len(kes & mapped) / len(kes), 1) if kes else 0.0})
    lv = pd.DataFrame(rows)
    # KE 1419 and KE 1420 carry no biological level: AOP-Wiki replaced them in AOP 232 with
    # KE 2421 and KE 2423 AFTER the organ networks were frozen, so they sit in the liver network
    # but not in the pinned release inventory. They are named rather than silently dropped,
    # because a level breakdown that quietly totals 536 against a stated 538 invites the reader
    # to think a row is missing.
    ORPHANS = {"KE 1419", "KE 1420"}
    missing = {k for k in union if level_of.get(k) is None}
    if missing != ORPHANS:
        raise SystemExit(f"Key Events with no biological level: {sorted(missing)}; "
                         f"expected exactly {sorted(ORPHANS)}")
    if lv["kes"].sum() != len(union) - len(ORPHANS):
        raise SystemExit(f"levels cover {lv['kes'].sum()} KEs, expected "
                         f"{len(union) - len(ORPHANS)}")
    lv.to_csv(DERIVED / "method3_coverage_by_biolevel.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 3.9),
                             gridspec_kw={"width_ratios": [1.05, 1], "wspace": 0.42})
    panel_a(axes[0], cov)
    panel_b(axes[1], lv)

    m1.source_note(fig, -0.13,
                   "Approved KE–WikiPathways mappings held by the molAOP Builder on 14 August "
                   "2026 (432 pairs over 314 Key Events, 158 distinct pathways), joined to the "
                   "four frozen organ networks. Counts include\nmappings that predate the D3.1 "
                   "and D3.2 curation campaigns. The four organ bars must not be summed: 50 of "
                   "the 538 Key Events belong to more than one network. Panel B covers 536 of "
                   "the 538 — KE 1419 and\nKE 1420 were replaced in AOP-Wiki after the networks "
                   "were frozen and carry no biological level in the pinned release.",
                   fontsize=7.6, linespacing=1.5)

    m1.save_figure(fig, "FIG10-method3-coverage")
    plt.close(fig)
    print(f"  -> {(DERIVED / 'method3_coverage_by_biolevel.tsv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
