#!/usr/bin/env python3
"""Does the proposal's cascade framing survive contact with the evidence? (Q7)

The full proposal describes the three methods as a **sequential** approach: text mining gives a
"first coarse gene annotation", "further improved for the underrepresented KEs" by the ontology
step, and finally molecular AOPs for "complex events either not captured by existing ontologies
or for which the previous step produced multiple matching".

That is a testable claim and the deliverable is obliged to answer it either way. It has two
parts, and they can come apart:

  A. **Does each step add coverage the previous ones lack?**  A waterfall of residuals.
  B. **Does each step preferentially rescue what the previous ones MISSED?**  This is the part
     that makes it a cascade rather than three independent samples. If Method 2 covers the Key
     Events Method 1 missed at the same rate as the ones Method 1 already had, then the methods
     are drawing from the same pool and the sequential framing has no support — even if A holds.

B is the discriminating test, and it is easy to get wrong. Coverage is confounded by curation
scope: Method 2 only covers Key Events a curator opened, and that choice was made by organ, not
by looking at Method 1's residual. So the conditional rates are computed BOTH on all 538 Key
Events and on the subset all three methods attempted, and only the latter is evidence about the
methods rather than about who curated what.

A note on what this can and cannot conclude. The three methods were in practice run in
PARALLEL, not in sequence — no step was ever pointed at the previous step's residual. So this
tests whether the data behave as the cascade predicts, not whether a cascade was executed. A
negative result means the framing does not describe these methods; it does not mean a genuine
cascade would fail.

Offline. Reads the comparison's own outputs.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method1"))
import m1lib as m1  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cmp_01_compare import METHODS, load  # noqa: E402

ROOT = m1.ROOT
OUT = ROOT / "data" / "comparison"
PROPOSAL_ORDER = ["method1", "method2", "method3"]
LABEL = {"method1": "Method 1 (text mining)",
         "method2": "Method 2 (GO BP)",
         "method3": "Method 3 (WikiPathways)"}


def two_by_two(a_hit: set, b_hit: set, universe: set) -> dict:
    """Coverage of b among the Key Events a covered, versus among the ones a missed."""
    a_yes, a_no = universe & a_hit, universe - a_hit
    r_yes = len(a_yes & b_hit) / len(a_yes) if a_yes else float("nan")
    r_no = len(a_no & b_hit) / len(a_no) if a_no else float("nan")
    out = {
        "n_previous_covered": len(a_yes),
        "n_previous_missed": len(a_no),
        "b_rate_where_previous_covered": round(100 * r_yes, 1),
        "b_rate_where_previous_missed": round(100 * r_no, 1),
        "rescue_difference_pp": round(100 * (r_no - r_yes), 1),
    }
    try:
        from scipy.stats import fisher_exact
        table = [[len(a_no & b_hit), len(a_no - b_hit)],
                 [len(a_yes & b_hit), len(a_yes - b_hit)]]
        odds, p = fisher_exact(table)
        out["odds_ratio_rescue_vs_not"] = round(float(odds), 3)
        out["fisher_p"] = float(f"{p:.3g}")
    except Exception as exc:  # scipy is not a hard dependency of this repo
        out["fisher_p"] = f"not computed ({exc})"
    return out


def waterfall(sets: dict, order: list[str], universe: set) -> list[dict]:
    rows, covered = [], set()
    for step, method in enumerate(order, 1):
        hit = {k for k in universe if sets[method].get(k)}
        new = hit - covered
        covered |= hit
        rows.append({
            "step": step,
            "method": method,
            "kes_covered_by_this_method": len(hit),
            "new_kes_added": len(new),
            "cumulative": len(covered),
            "cumulative_pct": round(100 * len(covered) / len(universe), 1),
            "residual_after": len(universe) - len(covered),
        })
    return rows


def main() -> None:
    m1.banner("Q7 — does the cascade framing survive the evidence?")

    sets, networks, attempted = load()
    universe = set().union(*networks.values())
    common = universe & attempted["method1"] & attempted["method2"] & attempted["method3"]
    hit = {m: {k for k in universe if sets[m].get(k)} for m in METHODS}

    # --- A. the waterfall, in the proposal's order --------------------------------------
    print(f"\n  A. Waterfall in the proposal's order, over {len(universe)} Key Events\n")
    wf = waterfall(sets, PROPOSAL_ORDER, universe)
    print(f"    {'step':<5} {'method':<10} {'covers':>7} {'NEW':>6} {'cumulative':>12} "
          f"{'residual':>9}")
    for r in wf:
        print(f"    {r['step']:<5} {r['method']:<10} {r['kes_covered_by_this_method']:>7} "
              f"{r['new_kes_added']:>6} {r['cumulative']:>7} ({r['cumulative_pct']:>4.1f} %) "
              f"{r['residual_after']:>9}")
    pd.DataFrame(wf).to_csv(OUT / "cascade_waterfall.tsv", sep="\t", index=False)

    # --- every ordering, to see whether the proposal's is the efficient one -------------
    orders = []
    for order in itertools.permutations(METHODS):
        rows = waterfall(sets, list(order), universe)
        orders.append({
            "order": " -> ".join(m[-1] for m in order),
            "is_proposal_order": list(order) == PROPOSAL_ORDER,
            "after_step_1": rows[0]["cumulative"],
            "after_step_2": rows[1]["cumulative"],
            "after_step_3": rows[2]["cumulative"],
            "step_2_marginal": rows[1]["new_kes_added"],
            "step_3_marginal": rows[2]["new_kes_added"],
        })
    od = pd.DataFrame(orders).sort_values("after_step_2", ascending=False)
    print(f"\n  Every ordering (the endpoint is identical by construction; what differs is how "
          f"fast it is reached)\n")
    print(f"    {'order':<12} {'step1':>6} {'step2':>6} {'step3':>6}   ")
    for _, r in od.iterrows():
        mark = "  <- proposal" if r["is_proposal_order"] else ""
        print(f"    {r['order']:<12} {r['after_step_1']:>6} {r['after_step_2']:>6} "
              f"{r['after_step_3']:>6}{mark}")
    od.to_csv(OUT / "cascade_orderings.tsv", sep="\t", index=False)

    # --- B. the discriminating test -----------------------------------------------------
    print("\n  B. Does each step preferentially rescue what the previous ones missed?")
    print("     (a cascade predicts the SECOND rate to be the higher one)\n")
    tests = []
    for universe_name, uni in (("all 538 KEs", universe),
                               ("KEs all three attempted", common)):
        for prev, nxt in (("method1", "method2"), ("method1", "method3")):
            t = two_by_two(hit[prev], hit[nxt], uni)
            t.update({"universe": universe_name, "previous": prev, "next": nxt})
            tests.append(t)
        # the real cascade step 3: everything the first two steps covered
        t = two_by_two(hit["method1"] | hit["method2"], hit["method3"], uni)
        t.update({"universe": universe_name, "previous": "method1 or method2",
                  "next": "method3"})
        tests.append(t)

    for t in tests:
        arrow = "rescues" if t["rescue_difference_pp"] > 0 else "does NOT rescue"
        print(f"    [{t['universe']:<24}] {t['next']} given {t['previous']}")
        print(f"        where previous covered: {t['b_rate_where_previous_covered']:>5.1f} %  "
              f"(n={t['n_previous_covered']})")
        print(f"        where previous missed : {t['b_rate_where_previous_missed']:>5.1f} %  "
              f"(n={t['n_previous_missed']})")
        print(f"        -> {arrow}, difference {t['rescue_difference_pp']:+.1f} pp, "
              f"p = {t.get('fisher_p')}")
    pd.DataFrame(tests).to_csv(OUT / "cascade_rescue_tests.tsv", sep="\t", index=False)

    # --- the residual: what no method reaches -------------------------------------------
    residual = universe - hit["method1"] - hit["method2"] - hit["method3"]
    inv = pd.read_csv(ROOT / "data/method1/raw/ke_inventory.csv")
    inv["keid"] = inv["keid"].str.strip()
    res = inv[inv["keid"].isin(residual)]
    allk = inv[inv["keid"].isin(universe)]
    # KE 1419 and KE 1420 sit in the frozen liver network but not in the pinned release:
    # AOP-Wiki replaced them in AOP 232 after the networks were frozen. They carry no biological
    # level, so the table below covers 536 of the 538. Named rather than silently dropped — a
    # breakdown that quietly totals 536 against a stated 538 reads as a missing row.
    orphans = sorted(universe - set(inv["keid"]))
    if orphans != ["KE 1419", "KE 1420"]:
        raise SystemExit(f"unexpected Key Events absent from the release inventory: {orphans}")
    print(f"\n  Residual: {len(residual)} of {len(universe)} Key Events are reached by no "
          f"method ({100 * len(residual) / len(universe):.1f} %)")
    print(f"    {len(orphans)} of the {len(universe)} carry no biological level in the pinned "
          f"release ({', '.join(orphans)}), so the table below covers "
          f"{len(universe) - len(orphans)}\n")
    lv = pd.DataFrame({
        "kes_in_networks": allk["biolevel"].value_counts(),
        "unreached": res["biolevel"].value_counts(),
    }).fillna(0).astype(int)
    lv["pct_unreached"] = (100 * lv["unreached"] / lv["kes_in_networks"]).round(1)
    lv = lv.sort_values("pct_unreached", ascending=False)
    print(f"    {'biological level':<16} {'in networks':>12} {'unreached':>10} {'%':>7}")
    for level, r in lv.iterrows():
        print(f"    {str(level):<16} {r['kes_in_networks']:>12} {r['unreached']:>10} "
              f"{r['pct_unreached']:>6.1f} %")
    lv.to_csv(OUT / "cascade_residual_by_level.tsv", sep="\t")
    res[["keid", "title", "biolevel"]].to_csv(
        OUT / "cascade_residual_kes.tsv", sep="\t", index=False)

    print("\n  -> data/comparison/cascade_waterfall.tsv")
    print("  -> data/comparison/cascade_orderings.tsv")
    print("  -> data/comparison/cascade_rescue_tests.tsv")
    print("  -> data/comparison/cascade_residual_by_level.tsv")
    print("  -> data/comparison/cascade_residual_kes.tsv")


if __name__ == "__main__":
    main()
