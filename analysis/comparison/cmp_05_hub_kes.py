#!/usr/bin/env python3
"""Are the Key Events no method reaches the central ones? (proposal risk table, T3.1)

The proposal's risk register carries one mitigation against "The AOP-Wiki does not fully capture
all sequences of KEs to lead to Adverse Outcomes": use complete organ-specific AOP networks, and
"we can prioritise the most central KEs (hub-KEs)".

D3.1 does the first half by construction — every result is computed on the four frozen organ
networks. The second half has never been answered, and it is answerable here: if the Key Events
the three methods miss are peripheral, the residual costs the project little; if they are hubs,
it costs a lot, because a hub Key Event carries the flow of several AOPs through one node.

Centrality is taken as degree in the union of the four networks' Key Event Relationships,
treated as undirected. Degree is the right measure for this question — a hub in the proposal's
sense is a Key Event many AOPs route through, which is exactly what a high KER count is.
Betweenness is reported alongside it when networkx is importable, as a second opinion.

Offline. Reads the frozen networks and the comparison's own coverage.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method1"))
import m1lib as m1  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cmp_01_compare import FOUR, METHODS, NODE_FILE, load  # noqa: E402

ROOT = m1.ROOT
NET = ROOT / "data" / "networks"
OUT = ROOT / "data" / "comparison"
EDGE_FILE = {"brain": "brain_edges.tsv", "kidney": "kidney_edges.tsv",
             "liver": "liver_edges.tsv", "lung": "lungs_edges.tsv"}


def union_edges() -> set[tuple[str, str]]:
    """Every KER appearing in any of the four networks, as an unordered pair."""
    edges: set[tuple[str, str]] = set()
    for organ in FOUR:
        df = pd.read_csv(NET / EDGE_FILE[organ], sep="\t")
        for src, tgt in zip(df["source"].str.strip(), df["target"].str.strip()):
            if src == tgt:
                continue
            edges.add((src, tgt) if src <= tgt else (tgt, src))
    return edges


def main() -> None:
    m1.banner("hub Key Events against what the three methods reach")

    sets, networks, attempted = load()
    union_kes = sorted(set().union(*networks.values()))
    edges = union_edges()

    degree = {ke: 0 for ke in union_kes}
    for a, b in edges:
        for node in (a, b):
            if node in degree:
                degree[node] += 1

    inventory = pd.read_csv(ROOT / "data" / "method1" / "raw" / "ke_inventory.csv")
    inventory["keid"] = inventory["keid"].astype(str).str.strip()
    title = dict(zip(inventory["keid"], inventory["title"].astype(str)))
    level = dict(zip(inventory["keid"], inventory["biolevel"].astype(str)))

    betweenness = {}
    try:
        import networkx as nx
        g = nx.Graph()
        g.add_nodes_from(union_kes)
        g.add_edges_from(edges)
        betweenness = nx.betweenness_centrality(g)
    except Exception as exc:
        print(f"  betweenness not computed ({exc})")

    rows = []
    for ke in union_kes:
        hits = [m for m in METHODS if sets[m].get(ke)]
        rows.append({
            "ke_id": ke,
            "title": title.get(ke, ""),
            "biolevel": level.get(ke, ""),
            "organs": ";".join(o for o in FOUR if ke in networks[o]),
            "degree": degree[ke],
            "betweenness": round(betweenness.get(ke, float("nan")), 6) if betweenness else "",
            "n_methods": len(hits),
            "covered": bool(hits),
        })
    df = pd.DataFrame(rows).sort_values(["degree", "ke_id"], ascending=[False, True])
    df.to_csv(OUT / "hub_kes.tsv", sep="\t", index=False)

    cov, res = df[df["covered"]], df[~df["covered"]]
    print(f"\n  {len(df)} Key Events, {len(edges)} distinct KERs in the four-network union")
    print(f"  reached by at least one method: {len(cov)}   reached by none: {len(res)}")
    print("\n  degree            median   mean   max   isolated (degree 0)")
    for name, part in (("reached", cov), ("unreached", res)):
        print(f"   {name:<16}{part['degree'].median():>6.1f}"
              f"{part['degree'].mean():>7.2f}{part['degree'].max():>6}"
              f"{int((part['degree'] == 0).sum()):>10}")

    # Coverage against centrality, in bands that a reader can act on.
    bands = [(0, 0, "isolated"), (1, 1, "degree 1"), (2, 3, "degree 2-3"),
             (4, 7, "degree 4-7"), (8, 10**6, "degree 8+")]
    band_rows = []
    for lo, hi, label in bands:
        part = df[(df["degree"] >= lo) & (df["degree"] <= hi)]
        if part.empty:
            continue
        band_rows.append({
            "band": label,
            "kes": len(part),
            "reached": int(part["covered"].sum()),
            "pct_reached": round(100 * part["covered"].mean(), 1),
            "mean_methods": round(part["n_methods"].mean(), 2),
        })
    bands_df = pd.DataFrame(band_rows)
    bands_df.to_csv(OUT / "hub_coverage_by_degree.tsv", sep="\t", index=False)
    print("\n" + bands_df.to_string(index=False))

    top = df.nlargest(int(round(0.1 * len(df))), "degree")
    print(f"\n  top decile by degree ({len(top)} KEs, degree >= {top['degree'].min()}): "
          f"{int(top['covered'].sum())} reached ({100 * top['covered'].mean():.1f} %), "
          f"{top['n_methods'].mean():.2f} methods each on average")
    print("\n  the ten most connected Key Events:")
    for _, r in df.head(10).iterrows():
        print(f"   deg {r['degree']:>3}  {r['n_methods']} method(s)  {r['ke_id']:<9} "
              f"{r['title'][:58]}")
    unreached_hubs = res[res["degree"] >= top["degree"].min()]
    print(f"\n  hub Key Events (top decile) that no method reaches: {len(unreached_hubs)}")
    for _, r in unreached_hubs.iterrows():
        print(f"   deg {r['degree']:>3}  {r['ke_id']:<9} {r['biolevel']:<12} {r['title'][:50]}")

    print(f"\n  wrote {OUT / 'hub_kes.tsv'}")
    print(f"  wrote {OUT / 'hub_coverage_by_degree.tsv'}")


if __name__ == "__main__":
    main()
