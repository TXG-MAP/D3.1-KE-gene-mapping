#!/usr/bin/env python3
"""How far the builder's gene sets are from the current upstream resources.

Methods 2 and 3 are reported from the molAOP Builder's own exports. Both were also derived
independently, straight from the upstream resources, and the gap between the two is a
reportable property of the deliverable rather than a discrepancy to resolve:

    Method 2   builder GO snapshot 2026-06-15   vs   GOA human 2026-05-21 + GO releases/2026-07-26
    Method 3   builder WP snapshot 2026-07-10   vs   WikiPathways RDF 20260810

The mappings are identical in both arms — the same curator approved the same KE-to-term and
KE-to-pathway pairs. What differs is only the gene membership of those terms and pathways, so
this measures corpus staleness alone, which is exactly what upstream molAOP-builder #263 is
open about.

Two reasons this belongs in the report rather than in a footnote. It bounds how much a reader
re-deriving these gene sets from WikiPathways or the GO today would get a different answer; and
it shows the direction of the effect, which is not symmetric between the two methods.

Offline. Reads four derived tables.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method1"))
import m1lib as m1  # noqa: E402

ROOT = m1.ROOT
OUT = ROOT / "data" / "comparison"

ARMS = {
    "method2": {
        "builder": ("data/method2/derived/method2_builder_ke_ensembl_sets.tsv", "ensembl", "|"),
        "upstream": ("data/method2/derived/method2_ke_gene_sets.tsv", "ensembl", "|"),
        "builder_pin": "builder GO snapshot 2026-06-15",
        "upstream_pin": "GOA human 2026-05-21 + GO releases/2026-07-26",
    },
    "method3": {
        "builder": ("data/method3/derived/method3_builder_ke_ensembl_sets.tsv", "ensembl", "|"),
        "upstream": ("data/method3/derived/method3_ke_ensembl_sets.tsv", "ensembl", "|"),
        "builder_pin": "builder WikiPathways snapshot 2026-07-10",
        "upstream_pin": "WikiPathways RDF 20260810",
    },
}


def load(spec) -> dict[str, set[str]]:
    path, col, sep = spec
    df = pd.read_csv(ROOT / path, sep="\t")
    out = {}
    for r in df.itertuples():
        value = getattr(r, col, "")
        out[r.ke_id.strip()] = ({v.strip() for v in value.split(sep) if v.strip()}
                                if isinstance(value, str) else set())
    return out


def main() -> None:
    m1.banner("builder vs upstream: how stale is the builder's corpus")

    rows, per_ke = [], []
    for method, spec in ARMS.items():
        b, u = load(spec["builder"]), load(spec["upstream"])
        shared = sorted(set(b) & set(u), key=lambda k: int(k.split()[-1]))
        both_nonempty = [k for k in shared if b[k] and u[k]]

        identical = sum(1 for k in both_nonempty if b[k] == u[k])
        b_genes = set().union(*b.values()) if b else set()
        u_genes = set().union(*u.values()) if u else set()
        sizes = [(len(b[k]), len(u[k])) for k in both_nonempty]

        rows.append({
            "method": method,
            "builder_pin": spec["builder_pin"],
            "upstream_pin": spec["upstream_pin"],
            "kes_in_both_arms": len(both_nonempty),
            "kes_with_identical_gene_sets": identical,
            "pct_identical": round(100 * identical / len(both_nonempty), 1)
            if both_nonempty else 0.0,
            "builder_genes": len(b_genes),
            "upstream_genes": len(u_genes),
            "genes_only_in_builder": len(b_genes - u_genes),
            "genes_only_upstream": len(u_genes - b_genes),
            "median_builder_set": float(pd.Series([s[0] for s in sizes]).median())
            if sizes else 0.0,
            "median_upstream_set": float(pd.Series([s[1] for s in sizes]).median())
            if sizes else 0.0,
            "median_pct_of_builder_set_still_upstream": round(100 * float(pd.Series(
                [len(b[k] & u[k]) / len(b[k]) for k in both_nonempty]).median()), 1)
            if both_nonempty else 0.0,
        })

        for k in both_nonempty:
            per_ke.append({
                "method": method, "ke_id": k,
                "builder_genes": len(b[k]), "upstream_genes": len(u[k]),
                "shared": len(b[k] & u[k]),
                "builder_only": len(b[k] - u[k]), "upstream_only": len(u[k] - b[k]),
            })

    drift = pd.DataFrame(rows)
    for _, r in drift.iterrows():
        print(f"\n  {r['method']}   {r['builder_pin']}")
        print(f"    vs           {r['upstream_pin']}")
        print(f"    {r['kes_in_both_arms']} KEs in both arms, "
              f"{r['kes_with_identical_gene_sets']} with identical gene sets "
              f"({r['pct_identical']} %)")
        print(f"    genes: {r['builder_genes']} builder / {r['upstream_genes']} upstream; "
              f"{r['genes_only_in_builder']} only in the builder, "
              f"{r['genes_only_upstream']} only upstream")
        print(f"    median set {r['median_builder_set']:.0f} vs "
              f"{r['median_upstream_set']:.0f} genes; "
              f"{r['median_pct_of_builder_set_still_upstream']} % of the median builder set is "
              f"still present upstream")

    drift.to_csv(OUT / "source_route_drift.tsv", sep="\t", index=False)
    pd.DataFrame(per_ke).to_csv(OUT / "source_route_drift_by_ke.tsv", sep="\t", index=False)
    print("\n  -> data/comparison/source_route_drift.tsv")
    print("  -> data/comparison/source_route_drift_by_ke.tsv")


if __name__ == "__main__":
    main()
