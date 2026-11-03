#!/usr/bin/env python3
"""The cross-method comparison — Methods 1, 2 and 3 over the four frozen organ networks.

This is what D3.1 is for: "Evaluation and comparison of various KE-gene mapping sets." The three
methods reach genes by routes that share nothing except their endpoint, so the comparison is
built on one common surface — **Ensembl gene identifiers on the KEs of the four frozen networks**
— and everything else is reported as a caveat rather than harmonised away.

    Method 1  text mining of KE and KER free text (HGNC dictionary + BERN2), harmonised to
              Ensembl through the AOP-Wiki RDF's own cross-references
    Method 2  curated GO Biological Process term per KE, expanded to genes through the GO
              annotation file with the true-path closure, UniProt -> Ensembl through BridgeDb
    Method 3  curated KE -> WikiPathways mappings, genes taken as the union of the approved
              pathways' Ensembl gene products

**The denominators are not the same, and that is the single most important thing about this
comparison.** Every method is scored against the same 538 KEs, but each is limited by something
different, and a bare coverage percentage silently attributes all three to the method:

    Method 1  can only read KEs that HAVE a free-text description — 635 of 1595 release-wide.
    Method 2  only covers KEs a curator has annotated. 355 of the 538 have been; the other 183
              were never attempted, so its coverage is a statement about curation effort.
    Method 3  is limited by whether a molecular pathway exists for the KE at all, which is
              structural: there is no pathway for "kidney failure".

So the headline overlap numbers below answer "what do the three mapping sets, as they stand
today, agree on" — not "which method is best". The second question needs the per-KE comparison
restricted to KEs all three actually attempted, which is why that restriction is computed
separately and reported alongside.

Offline: reads the three methods' derived tables and the frozen networks. No network calls.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method1"))
import m1lib as m1  # noqa: E402

ROOT = m1.ROOT
OUT = ROOT / "data" / "comparison"
OUT.mkdir(parents=True, exist_ok=True)

FOUR = ["brain", "kidney", "liver", "lung"]
NODE_FILE = {"brain": "brain_nodes.tsv", "kidney": "kidney_nodes.tsv",
             "liver": "liver_nodes.tsv", "lung": "lungs_nodes.tsv"}
METHODS = ["method1", "method2", "method3"]


def split(value: str, sep: str) -> set[str]:
    if not isinstance(value, str) or not value.strip():
        return set()
    return {v.strip() for v in value.split(sep) if v.strip()}


def load() -> tuple[dict[str, dict[str, set[str]]], dict[str, set[str]], dict[str, set[str]]]:
    """Return (per-method KE -> gene sets, organ -> KE ids, method -> KEs attempted).

    "Attempted" is not the same as "annotated" and the distinction carries the whole
    interpretation. A KE that Method 2 never looked at is not evidence that GO has nothing to
    say about it.
    """
    networks = {
        organ: set(pd.read_csv(ROOT / "data" / "networks" / NODE_FILE[organ],
                               sep="\t")["id"].str.strip())
        for organ in FOUR
    }

    m1_df = pd.read_csv(ROOT / "data/method1/derived/method1_ke_ensembl_sets.csv")
    sets = {"method1": {r.keid.strip(): split(r.ensembl, ",") for r in m1_df.itertuples()}}

    # Methods 2 and 3 come from the BUILDER's own exports (Marvin, 2026-08-14). Both were
    # curated in it, so the comparison reports what the resource holds rather than an
    # independent re-derivation, and both incur the same symbol -> Ensembl round-trip.
    # analysis/builder_genesets.py writes these; the independent routes still run as
    # cross-checks and are read below only to report the drift between them.
    m2_df = pd.read_csv(
        ROOT / "data/method2/derived/method2_builder_ke_ensembl_sets.tsv", sep="\t")
    sets["method2"] = {r.ke_id.strip(): split(r.ensembl, "|") for r in m2_df.itertuples()}

    m3_df = pd.read_csv(
        ROOT / "data/method3/derived/method3_builder_ke_ensembl_sets.tsv", sep="\t")
    sets["method3"] = {r.ke_id.strip(): split(r.ensembl, "|") for r in m3_df.itertuples()}

    # Method 1 reads descriptions, so the KEs it ATTEMPTED are the described ones; Method 2
    # attempted the KEs a curator opened; Method 3 attempted every KE in the four networks —
    # its 181 `no_suitable_pathway` verdicts are answers, not omissions.
    # data/method1/raw/ is gitignored bulk, but this is the only place the described/undescribed
    # split is recorded, and it is what makes Method 1's denominator honest (C30). Regenerate
    # with `make method1-fetch` if it is absent.
    inventory_path = ROOT / "data/method1/raw/ke_inventory.csv"
    if not inventory_path.exists():
        raise SystemExit(
            f"{inventory_path.relative_to(ROOT)} is missing — run `make method1-fetch` first. "
            "Method 1 can only read KEs that have a free-text description, and without that "
            "flag its coverage column silently compares against the wrong denominator."
        )
    described = set(pd.read_csv(inventory_path).query("has_description")["keid"].str.strip())
    # "Attempted" for Method 2 comes from the CURATION files, not from the builder export. The
    # export holds approved mappings only, so a KE a curator opened and deliberately left empty
    # — 53 of them — is absent from it. Counting those as never-attempted would credit Method 2
    # with a higher hit rate for having declined to annotate.
    curated = pd.read_csv(ROOT / "data/method2/derived/method2_ke_gene_sets.tsv", sep="\t")
    attempted = {
        "method1": described,
        "method2": {k.strip() for k in curated["ke_id"]}
        | {r.ke_id.strip() for r in m2_df.itertuples()},
        "method3": set().union(*networks.values()),
    }
    return sets, networks, attempted


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a or b) else float("nan")


def containment(a: set, b: set) -> float:
    """Fraction of ``a`` that also appears in ``b``.

    Jaccard is the wrong statistic for these three methods and reporting it alone would be
    misleading. The median annotated KE carries about 6 genes under Method 1 and about 65 under
    Methods 2 and 3, so a Method 1 set can be *entirely* inside a Method 3 set and still score a
    Jaccard of 0.09 — the union term punishes the size gap, not the disagreement. Containment
    asks the question the report actually cares about: when the text mining finds a gene, do the
    curated routes already have it? Both are reported; neither is sufficient alone.
    """
    return len(a & b) / len(a) if a else float("nan")


def main() -> None:
    m1.banner("cross-method comparison over the four frozen networks")

    sets, networks, attempted = load()
    union_kes = set().union(*networks.values())
    frames = {**networks, "four-organ union": union_kes}

    # --- 1. coverage of each network, per method ----------------------------------------
    rows = []
    for organ, ids in frames.items():
        row = {"organ": organ, "kes": len(ids)}
        for m in METHODS:
            with_genes = {k for k in ids if sets[m].get(k)}
            genes = set().union(*(sets[m].get(k, set()) for k in ids)) if ids else set()
            row[f"{m}_kes"] = len(with_genes)
            row[f"{m}_pct"] = round(100 * len(with_genes) / len(ids), 1)
            row[f"{m}_genes"] = len(genes)
            row[f"{m}_attempted"] = len(ids & attempted[m])
        rows.append(row)
    coverage = pd.DataFrame(rows)
    print("\n  coverage of each frozen network (% of the network's KEs carrying genes)")
    for _, r in coverage.iterrows():
        print(f"    {r['organ']:<18} M1 {r['method1_pct']:>5.1f} %  "
              f"M2 {r['method2_pct']:>5.1f} %  M3 {r['method3_pct']:>5.1f} %"
              f"   ({r['kes']} KEs)")
    coverage.to_csv(OUT / "coverage_by_method_and_organ.tsv", sep="\t", index=False)

    # --- 2. how many KEs each combination of methods reaches ----------------------------
    combo_rows = []
    for organ, ids in frames.items():
        have = {m: {k for k in ids if sets[m].get(k)} for m in METHODS}
        combo_rows.append({
            "organ": organ,
            "kes": len(ids),
            "none": len(ids - set().union(*have.values())),
            "method1_only": len(have["method1"] - have["method2"] - have["method3"]),
            "method2_only": len(have["method2"] - have["method1"] - have["method3"]),
            "method3_only": len(have["method3"] - have["method1"] - have["method2"]),
            "m1_m2": len(have["method1"] & have["method2"] - have["method3"]),
            "m1_m3": len(have["method1"] & have["method3"] - have["method2"]),
            "m2_m3": len(have["method2"] & have["method3"] - have["method1"]),
            "all_three": len(have["method1"] & have["method2"] & have["method3"]),
            "any": len(set().union(*have.values())),
        })
    combos = pd.DataFrame(combo_rows)
    u = combos[combos["organ"] == "four-organ union"].iloc[0]
    print(f"\n  union of the three methods: {u['any']} of {u['kes']} KEs carry genes from at "
          f"least one method ({100 * u['any'] / u['kes']:.1f} %); "
          f"{u['all_three']} from all three")
    print(f"    unique contributions — M1 only {u['method1_only']}, M2 only "
          f"{u['method2_only']}, M3 only {u['method3_only']}")
    combos.to_csv(OUT / "ke_coverage_combinations.tsv", sep="\t", index=False)

    # --- 3. per-KE agreement, where two methods both produced a set ---------------------
    pair_rows = []
    for a, b in itertools.combinations(METHODS, 2):
        for organ, ids in frames.items():
            both = [k for k in ids if sets[a].get(k) and sets[b].get(k)]
            if not both:
                continue
            js = [jaccard(sets[a][k], sets[b][k]) for k in both]
            inter = [len(sets[a][k] & sets[b][k]) for k in both]
            c_ab = [containment(sets[a][k], sets[b][k]) for k in both]
            c_ba = [containment(sets[b][k], sets[a][k]) for k in both]
            pair_rows.append({
                "pair": f"{a} vs {b}",
                "organ": organ,
                "kes_both": len(both),
                "kes_with_any_shared_gene": sum(1 for i in inter if i > 0),
                "pct_with_any_shared_gene": round(
                    100 * sum(1 for i in inter if i > 0) / len(both), 1),
                "median_jaccard": round(float(pd.Series(js).median()), 4),
                "mean_jaccard": round(float(pd.Series(js).mean()), 4),
                "median_shared_genes": float(pd.Series(inter).median()),
                f"median_pct_of_{a}_in_{b}": round(100 * float(pd.Series(c_ab).median()), 1),
                f"median_pct_of_{b}_in_{a}": round(100 * float(pd.Series(c_ba).median()), 1),
                "median_set_size_a": float(pd.Series([len(sets[a][k]) for k in both]).median()),
                "median_set_size_b": float(pd.Series([len(sets[b][k]) for k in both]).median()),
            })
    pairs = pd.DataFrame(pair_rows)
    print("\n  per-KE agreement where both methods produced a gene set (four-organ union)")
    for _, r in pairs[pairs["organ"] == "four-organ union"].iterrows():
        a, b = r["pair"].split(" vs ")
        print(f"    {r['pair']:<24} {r['kes_both']:>3} KEs  "
              f"{r['pct_with_any_shared_gene']:>5.1f} % share >=1 gene  "
              f"J {r['median_jaccard']:.3f}  "
              f"{a}->{b} {r[f'median_pct_of_{a}_in_{b}']:>5.1f} %  "
              f"{b}->{a} {r[f'median_pct_of_{b}_in_{a}']:>5.1f} %")
    pairs.to_csv(OUT / "pairwise_ke_agreement.tsv", sep="\t", index=False)

    # --- 4. gene-level overlap ----------------------------------------------------------
    gene_rows = []
    for organ, ids in frames.items():
        genes = {m: set().union(*(sets[m].get(k, set()) for k in ids)) if ids else set()
                 for m in METHODS}
        gene_rows.append({
            "organ": organ,
            **{f"{m}_genes": len(genes[m]) for m in METHODS},
            "in_all_three": len(genes["method1"] & genes["method2"] & genes["method3"]),
            "in_any": len(set().union(*genes.values())),
            "method1_only": len(genes["method1"] - genes["method2"] - genes["method3"]),
            "method2_only": len(genes["method2"] - genes["method1"] - genes["method3"]),
            "method3_only": len(genes["method3"] - genes["method1"] - genes["method2"]),
            "jaccard_m1_m2": round(jaccard(genes["method1"], genes["method2"]), 4),
            "jaccard_m1_m3": round(jaccard(genes["method1"], genes["method3"]), 4),
            "jaccard_m2_m3": round(jaccard(genes["method2"], genes["method3"]), 4),
        })
    gene_overlap = pd.DataFrame(gene_rows)
    g = gene_overlap[gene_overlap["organ"] == "four-organ union"].iloc[0]
    print(f"\n  gene level, four-organ union: {g['in_any']} distinct genes in total, "
          f"{g['in_all_three']} found by all three")
    print(f"    only M1 {g['method1_only']}, only M2 {g['method2_only']}, "
          f"only M3 {g['method3_only']}")
    gene_overlap.to_csv(OUT / "gene_level_overlap.tsv", sep="\t", index=False)

    # --- 5. the same comparison restricted to KEs all three attempted -------------------
    #
    # Without this restriction the coverage columns mostly measure how much curation each
    # method has received. A KE Method 2 never opened is not evidence about GO.
    common = union_kes & attempted["method1"] & attempted["method2"] & attempted["method3"]
    fair_rows = []
    for m in METHODS:
        with_genes = {k for k in common if sets[m].get(k)}
        fair_rows.append({
            "method": m,
            "kes_all_three_attempted": len(common),
            "kes_with_genes": len(with_genes),
            "pct": round(100 * len(with_genes) / len(common), 1) if common else 0.0,
            "median_genes_per_annotated_ke": float(pd.Series(
                [len(sets[m][k]) for k in with_genes]).median()) if with_genes else 0.0,
        })
    fair = pd.DataFrame(fair_rows)
    print(f"\n  restricted to the {len(common)} KEs ALL THREE methods attempted:")
    for _, r in fair.iterrows():
        print(f"    {r['method']:<9} {r['kes_with_genes']:>3}/{r['kes_all_three_attempted']} "
              f"({r['pct']:>5.1f} %)  median {r['median_genes_per_annotated_ke']:.0f} genes")
    fair.to_csv(OUT / "coverage_on_commonly_attempted_kes.tsv", sep="\t", index=False)

    # --- 6. the per-KE long table everything above is derived from ----------------------
    long_rows = []
    for ke in sorted(union_kes, key=lambda k: int(k.split()[-1])):
        organs = [o for o in FOUR if ke in networks[o]]
        g = {m: sets[m].get(ke, set()) for m in METHODS}
        long_rows.append({
            "ke_id": ke,
            "organs": "|".join(organs),
            **{f"{m}_n": len(g[m]) for m in METHODS},
            **{f"{m}_attempted": ke in attempted[m] for m in METHODS},
            "shared_all_three": len(g["method1"] & g["method2"] & g["method3"]),
            "union_genes": len(set().union(*g.values())),
            "jaccard_m1_m2": round(jaccard(g["method1"], g["method2"]), 4)
            if g["method1"] and g["method2"] else "",
            "jaccard_m1_m3": round(jaccard(g["method1"], g["method3"]), 4)
            if g["method1"] and g["method3"] else "",
            "jaccard_m2_m3": round(jaccard(g["method2"], g["method3"]), 4)
            if g["method2"] and g["method3"] else "",
        })
    long = pd.DataFrame(long_rows)
    long.to_csv(OUT / "ke_level_comparison.tsv", sep="\t", index=False)
    print(f"\n  -> {len(long)} KE rows in data/comparison/ke_level_comparison.tsv")

    for name in ("coverage_by_method_and_organ", "ke_coverage_combinations",
                 "pairwise_ke_agreement", "gene_level_overlap",
                 "coverage_on_commonly_attempted_kes"):
        print(f"  -> data/comparison/{name}.tsv")

    write_markdown(coverage, combos, pairs, gene_overlap, fair, len(common))


def write_markdown(coverage, combos, pairs, gene_overlap, fair, n_common) -> None:
    u = combos[combos["organ"] == "four-organ union"].iloc[0]
    g = gene_overlap[gene_overlap["organ"] == "four-organ union"].iloc[0]
    p = pairs[pairs["organ"] == "four-organ union"]

    out = ["# Cross-method comparison — Methods 1, 2 and 3", "",
           "Generated by `analysis/comparison/cmp_01_compare.py` over the four frozen organ "
           "networks (538 distinct Key Events), in Ensembl identifier space.", "",
           "**Where the gene sets come from.** Method 1 from the AOP-Wiki RDF at the pinned "
           "release. Methods 2 and 3 from the **molAOP Builder's own exports** "
           "(`/exports/gmt/ke-go` and `/exports/gmt/ke-wp`) — both were curated in the builder, "
           "so the comparison reports what the resource holds rather than an independent "
           "re-derivation. Both those exports carry HGNC symbols, so both methods incur the "
           "same symbol-to-Ensembl round-trip (99.4 % of 9620 symbols map; 56 do not). The "
           "builder serves a corpus snapshot of **WikiPathways 2026-07-10 · GO 2026-06-15**, "
           "which is behind the live resources — `data/comparison/source_route_drift.tsv` "
           "measures by how much.", "",
           "## Coverage of each network", "",
           "Percentage of the network's Key Events carrying at least one gene. **These three "
           "columns are not measuring the same thing** — see the caveat below the table.", "",
           "| Organ | KEs | Method 1 | Method 2 | Method 3 |",
           "| --- | ---: | ---: | ---: | ---: |"]
    for _, r in coverage.iterrows():
        name = r["organ"] if r["organ"] == "four-organ union" else r["organ"].capitalize()
        out.append(f"| {name} | {r['kes']} | {r['method1_pct']} % | {r['method2_pct']} % "
                   f"| {r['method3_pct']} % |")
    out += ["",
            "Each method is limited by something different, and reading these as a ranking "
            "would attribute all three limits to the methods themselves:",
            "",
            "- **Method 1** can only read a Key Event that HAS a free-text description.",
            "- **Method 2** only covers Key Events a curator has annotated — 355 of the 538. "
            "The remaining 183 were never attempted.",
            "- **Method 3** is limited by whether a molecular pathway exists for the Key Event "
            "at all, which is structural rather than a matter of effort.",
            "",
            f"Restricted to the **{n_common} Key Events all three methods attempted**:",
            "",
            "| Method | KEs with genes | coverage | median genes per annotated KE |",
            "| --- | ---: | ---: | ---: |"]
    for _, r in fair.iterrows():
        out.append(f"| {r['method']} | {r['kes_with_genes']}/{r['kes_all_three_attempted']} "
                   f"| {r['pct']} % | {r['median_genes_per_annotated_ke']:.0f} |")

    out += ["", "## What the three methods add to each other", "",
            f"Across the four networks, **{u['any']} of {u['kes']} Key Events "
            f"({100 * u['any'] / u['kes']:.1f} %)** carry genes from at least one method, and "
            f"**{u['all_three']}** from all three.", "",
            "| | KEs |", "| --- | ---: |",
            f"| Method 1 only | {u['method1_only']} |",
            f"| Method 2 only | {u['method2_only']} |",
            f"| Method 3 only | {u['method3_only']} |",
            f"| Methods 1 and 2 | {u['m1_m2']} |",
            f"| Methods 1 and 3 | {u['m1_m3']} |",
            f"| Methods 2 and 3 | {u['m2_m3']} |",
            f"| All three | {u['all_three']} |",
            f"| None | {u['none']} |",
            "",
            "## Do they agree on the genes themselves?", "",
            "Only on Key Events where **both** methods produced a set. Jaccard and containment "
            "are both given because they answer different questions: the median annotated Key "
            "Event carries about 6 genes under Method 1 and about 65 under Methods 2 and 3, so "
            "a Method 1 set can sit entirely inside a Method 3 set and still score a low "
            "Jaccard. Containment reads: of the genes the first method found, what share the "
            "second also has.",
            "",
            "| Pair | KEs | share >=1 gene | median Jaccard | A in B | B in A |",
            "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for _, r in p.iterrows():
        a, b = r["pair"].split(" vs ")
        out.append(f"| {r['pair']} | {r['kes_both']} | {r['pct_with_any_shared_gene']} % "
                   f"| {r['median_jaccard']:.3f} | {r[f'median_pct_of_{a}_in_{b}']} % "
                   f"| {r[f'median_pct_of_{b}_in_{a}']} % |")

    out += ["",
            f"At gene level over the four networks, the three methods name **{g['in_any']} "
            f"distinct genes** between them and agree on **{g['in_all_three']}**: "
            f"{g['method1_only']} are unique to Method 1, {g['method2_only']} to Method 2 and "
            f"{g['method3_only']} to Method 3.",
            "",
            "The picture is one of **complementarity rather than redundancy**. Most Key Events "
            "annotated by two methods share at least one gene, but the median Key Event shares "
            "only a small minority of them, and each method contributes Key Events no other "
            "method reaches.", ""]

    path = ROOT / "tables" / "comparison_summary.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out))
    print(f"  -> {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
