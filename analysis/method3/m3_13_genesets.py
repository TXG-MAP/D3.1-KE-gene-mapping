#!/usr/bin/env python3
"""Stage 13 — turn the approved KE->WikiPathways mappings into KE->gene sets in Ensembl.

Method 3 does not produce genes. It produces *pathways*, and a KE's gene set is the union of
the gene sets of the pathways approved for it. This stage performs that union, which is what
makes Method 3 comparable with Methods 1 and 2 at all.

Genes come from **WikiPathways directly** (SPARQL, release pinned), not from the builder's GMT
export, for two reasons set out in queries/method3a_pathway_genes.rq: the GMT carries ambiguous
HGNC symbols where WikiPathways holds the Ensembl cross-reference outright, and the builder's
corpus is an unpinned snapshot. The GMT is still compared against, and the drift is written out
— it is a reportable fact about Method 3's reproducibility.

**The circularity caveat travels with every output of this stage.** Method 3's gene column is
not independent evidence for a KE-pathway mapping: the pathway was chosen partly on its genes,
so its genes cannot then corroborate it (Method 3 limitation 4; the Ensembl sets are 92-100 %
contained in the KE's own mapped pathways by construction). These sets are legitimate as *the
output of the method* — which is what the cross-method comparison consumes — and illegitimate
as evidence *for* the method.

Network: WikiPathways SPARQL in batches, plus one GET to the builder for the GMT comparison.
"""
from __future__ import annotations

import argparse
import json
import re
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import date

import pandas as pd

import m3lib as L

WP_SPARQL = "https://sparql.wikipathways.org/sparql"

#: The WikiPathways RDF release this stage is pinned to. Read from the endpoint's own void
#: metadata and asserted before any gene row is used — the endpoint reloads on its own schedule
#: and a silent move would change every gene set below without changing a single line of code.
WP_RELEASE = "WikiPathways RDF 20260810"
WP_CREATED_ON = "2026-08-10"

QUERY = L.ROOT / "queries" / "method3a_pathway_genes.rq"
BATCH = 40
STAMP = date.today().isoformat()

FOUR = ["brain", "kidney", "liver", "lung"]
NODE_FILE = {"brain": "brain_nodes.tsv", "kidney": "kidney_nodes.tsv",
             "liver": "liver_nodes.tsv", "lung": "lungs_nodes.tsv"}


def sparql(query: str, timeout: int = 180) -> pd.DataFrame:
    url = WP_SPARQL + "?" + urllib.parse.urlencode({"query": query})
    req = urllib.request.Request(url, headers={
        "Accept": "application/sparql-results+json",
        "User-Agent": L.USER_AGENT,
    })
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.loads(resp.read())
    cols = payload["head"]["vars"]
    return pd.DataFrame([
        {c: b.get(c, {}).get("value") for c in cols}
        for b in payload["results"]["bindings"]
    ], columns=cols)


def verify_wp_release() -> dict:
    """Refuse to run against a WikiPathways release other than the pinned one.

    Same failure mode as the AOP-Wiki endpoint, and the same treatment: the data behind the
    number must be identifiable months later. --allow-drift is available for a deliberate
    re-pin and prints loudly.
    """
    df = sparql("""
        PREFIX void: <http://rdfs.org/ns/void#>
        PREFIX pav:  <http://purl.org/pav/>
        PREFIX dcterms: <http://purl.org/dc/terms/>
        SELECT ?title ?created WHERE {
          ?d a void:Dataset ; dcterms:title ?title ; pav:createdOn ?created .
          FILTER(CONTAINS(STR(?title), "WikiPathways RDF"))
        }
    """)
    if df.empty:
        raise SystemExit("WikiPathways endpoint returned no dataset metadata")
    title = df["title"].iloc[0]
    created = str(df["created"].iloc[0])[:10]
    print(f"  WikiPathways release: {title} (created {created})")
    return {"title": title, "created_on": created,
            "matches_pin": title == WP_RELEASE and created == WP_CREATED_ON}


def fetch_pathway_genes(wpids: list[str]) -> pd.DataFrame:
    template = QUERY.read_text()
    frames = []
    for start in range(0, len(wpids), BATCH):
        block = wpids[start:start + BATCH]
        values = " ".join(f'"{w}"' for w in block)
        # The committed query carries a __WPIDS__ placeholder so the file on disk is the query
        # that ran, modulo the VALUES block. Appendix A reproduces it as written.
        frames.append(sparql(template.replace("__WPIDS__", values)))
        print(f"    {min(start + BATCH, len(wpids)):4d}/{len(wpids)} pathways", end="\r")
    print(" " * 40, end="\r")
    df = pd.concat(frames, ignore_index=True)
    df["ensembl"] = df["ensembl"].str.rsplit("/", n=1).str[-1]
    return df


def gmt_symbols() -> dict[str, set[str]]:
    """The builder's own gene sets per pathway, as symbols — for the drift comparison.

    ``.strip()`` on every field is not cosmetic. The last gene of each GMT line keeps its line
    terminator, so a naive split makes ``ADIPOR1\\n`` a symbol that matches nothing, and every
    pathway then reports one spurious builder-only gene.
    """
    body, _ = L.get("/exports/gmt/ke-wp", raw=True, pace=1.0)
    out: dict[str, set[str]] = defaultdict(set)
    for line in body.decode().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split("\t")]
        m = re.match(r"^KE(\d+)_.*_(WP\d+)$", parts[0])
        if m:
            out[m.group(2)] |= {g for g in parts[2:] if g}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-drift", action="store_true",
                    help="proceed even if the WikiPathways release has moved off the pin")
    args = ap.parse_args()

    L.banner("stage 13: Method 3 KE gene sets, via WikiPathways")

    release = verify_wp_release()
    if not release["matches_pin"]:
        msg = (f"WikiPathways has moved off the pin: expected '{WP_RELEASE}' created "
               f"{WP_CREATED_ON}, endpoint serves '{release['title']}' created "
               f"{release['created_on']}")
        if not args.allow_drift:
            raise SystemExit(
                msg + "\n\nRe-pin WP_RELEASE/WP_CREATED_ON deliberately and re-run every "
                      "Method 3 gene number, or pass --allow-drift and expect the numbers not "
                      "to match anything already reported."
            )
        print(f"  ! {msg}\n  ! continuing under --allow-drift")

    pairs = pd.read_csv(L.DERIVED_DIR / "method3_approved_pairs.tsv", sep="\t")
    frames = {o: pd.read_csv(L.ROOT / "data" / "networks" / NODE_FILE[o], sep="\t")
              for o in FOUR}
    union_kes = set().union(*(set(df["id"].str.strip()) for df in frames.values()))
    pairs = pairs[pairs["ke_id"].isin(union_kes)]
    wpids = sorted(pairs["wp_id"].unique(), key=lambda w: int(w[2:]))
    print(f"  {len(pairs)} approved pairs on the four networks, "
          f"{len(wpids)} distinct pathways")

    print("\n  fetching pathway gene sets")
    genes = fetch_pathway_genes(wpids)
    per_pathway = genes.groupby("wpid")["ensembl"].nunique()
    missing = [w for w in wpids if w not in set(genes["wpid"])]
    print(f"    {genes['ensembl'].nunique()} distinct Ensembl genes over "
          f"{per_pathway.size} pathways "
          f"(median {per_pathway.median():.0f} genes/pathway)")
    if missing:
        # A pathway with no human gene product is not an error: WikiPathways holds metabolic
        # and signalling maps whose nodes are metabolites. It does mean the KEs mapped only to
        # such pathways get an EMPTY gene set, which the comparison must treat as "Method 3
        # made a mapping but contributes no genes" rather than as "no mapping".
        print(f"    ! {len(missing)} pathways returned no human gene product: "
              f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}")

    L.save_table(genes.sort_values(["wpid", "ensembl"]),
                 L.DERIVED_DIR / "method3_pathway_genes.tsv")

    # --- KE gene sets -------------------------------------------------------------------
    by_pathway = genes.groupby("wpid")["ensembl"].apply(set).to_dict()
    rows = []
    for ke, grp in pairs.groupby("ke_id"):
        ensembl: set[str] = set()
        for wp in grp["wp_id"]:
            ensembl |= by_pathway.get(wp, set())
        rows.append({
            "ke_id": ke,
            "n_pathways": grp["wp_id"].nunique(),
            "n_genes": len(ensembl),
            "ensembl": "|".join(sorted(ensembl)),
        })
    ke_sets = pd.DataFrame(rows).sort_values("ke_id", key=lambda s: s.str.extract(r"(\d+)")[0].astype(int))
    empty = int((ke_sets["n_genes"] == 0).sum())
    print(f"\n  {len(ke_sets)} KEs carry a gene set, {ke_sets['n_genes'].sum()} KE-gene pairs, "
          f"{len(set('|'.join(ke_sets['ensembl']).split('|')) - {''})} distinct genes")
    if empty:
        print(f"    {empty} of them are EMPTY — mapped to pathways with no human gene product")
    L.save_table(ke_sets, L.DERIVED_DIR / "method3_ke_ensembl_sets.tsv")

    # --- per organ ----------------------------------------------------------------------
    organ_rows = []
    for organ in FOUR:
        ids = set(frames[organ]["id"].str.strip())
        sub = ke_sets[ke_sets["ke_id"].isin(ids)]
        gene_union = set("|".join(sub["ensembl"]).split("|")) - {""}
        organ_rows.append({
            "organ": organ,
            "kes_in_frame": len(ids),
            "kes_with_genes": int((sub["n_genes"] > 0).sum()),
            "pct_kes_with_genes": round(100 * (sub["n_genes"] > 0).sum() / len(ids), 1),
            "unique_genes": len(gene_union),
            "ke_gene_pairs": int(sub["n_genes"].sum()),
            "median_genes_per_annotated_ke": float(sub.loc[sub["n_genes"] > 0, "n_genes"].median())
            if (sub["n_genes"] > 0).any() else 0.0,
        })
    all_genes = set("|".join(ke_sets["ensembl"]).split("|")) - {""}
    organ_rows.append({
        "organ": "four-organ union",
        "kes_in_frame": len(union_kes),
        "kes_with_genes": int((ke_sets["n_genes"] > 0).sum()),
        "pct_kes_with_genes": round(100 * (ke_sets["n_genes"] > 0).sum() / len(union_kes), 1),
        "unique_genes": len(all_genes),
        "ke_gene_pairs": int(ke_sets["n_genes"].sum()),
        "median_genes_per_annotated_ke": float(
            ke_sets.loc[ke_sets["n_genes"] > 0, "n_genes"].median()),
    })
    organs = pd.DataFrame(organ_rows)
    print()
    for _, r in organs.iterrows():
        print(f"    {r['organ']:<18} {r['kes_with_genes']:>4}/{r['kes_in_frame']:<4} KEs "
              f"({r['pct_kes_with_genes']:>5.1f} %)  {r['unique_genes']:>5} genes")
    L.save_table(organs, L.DERIVED_DIR / "method3_genes_by_organ.tsv")

    # --- drift against the builder's own export ------------------------------------------
    # --- drift against the builder's own export ------------------------------------------
    #
    # This comparison is easy to get wrong in a way that manufactures an alarming number, and
    # the first version of this stage did. Two asymmetries have to be handled:
    #
    #   * The builder exports one SYMBOL per gene. WikiPathways' rdfs:label is not unique per
    #     Ensembl identifier — 2413 of 8430 gene products here carry more than one label
    #     (ENSG00000129965 is both "INS" and "Insulin"). Comparing set sizes therefore compares
    #     labels against genes and reports drift on almost every pathway. Sizes are compared on
    #     DISTINCT ENSEMBL; labels are used only to ask whether a builder symbol is findable.
    #   * A builder symbol absent from every WikiPathways label for that pathway is the signal
    #     that actually matters — it means the builder's snapshot holds a gene the current
    #     release does not, i.e. the corpus is stale (upstream molAOP-builder #263).
    print("\n  comparing against the builder's GMT")
    builder = gmt_symbols()
    wp_labels = genes.dropna(subset=["symbol"]).groupby("wpid")["symbol"].apply(set).to_dict()
    wp_ensembl = genes.groupby("wpid")["ensembl"].nunique().to_dict()
    drift = []
    for wp in wpids:
        b = builder.get(wp, set())
        unfindable = b - wp_labels.get(wp, set())
        drift.append({
            "wp_id": wp,
            "builder_symbols": len(b),
            "wikipathways_genes": wp_ensembl.get(wp, 0),
            "size_delta": wp_ensembl.get(wp, 0) - len(b),
            "builder_symbols_not_in_wikipathways": len(unfindable),
            "unfindable": "|".join(sorted(unfindable)),
        })
    drift_df = pd.DataFrame(drift)
    same_size = int((drift_df["size_delta"] == 0).sum())
    with_unfindable = int((drift_df["builder_symbols_not_in_wikipathways"] > 0).sum())
    print(f"    {same_size} of {len(wpids)} pathways hold the same number of genes in both; "
          f"median delta {drift_df['size_delta'].median():+.0f}")
    print(f"    {with_unfindable} pathways carry at least one builder symbol absent from "
          f"WikiPathways {WP_CREATED_ON} "
          f"({int(drift_df['builder_symbols_not_in_wikipathways'].sum())} symbols in total)")
    L.save_table(drift_df.sort_values("size_delta"),
                 L.DERIVED_DIR / "method3_corpus_drift.tsv")

    L.write_manifest("13_genesets", {
        "purpose": "KE -> Ensembl gene sets for Method 3, from the approved pathway mappings",
        "wikipathways_release": release,
        "wikipathways_pin": {"title": WP_RELEASE, "created_on": WP_CREATED_ON},
        "wikipathways_endpoint": WP_SPARQL,
        "query": "queries/method3a_pathway_genes.rq",
        "approved_pairs_used": len(pairs),
        "distinct_pathways": len(wpids),
        "pathways_without_human_gene_products": missing,
        "distinct_ensembl_genes": len(all_genes),
        "kes_with_a_gene_set": int((ke_sets["n_genes"] > 0).sum()),
        "kes_mapped_but_empty": empty,
        "ke_gene_pairs": int(ke_sets["n_genes"].sum()),
        "by_organ": organs.to_dict("records"),
        "corpus_drift": {
            "pathways_same_gene_count": same_size,
            "median_size_delta": float(drift_df["size_delta"].median()),
            "pathways_with_unfindable_builder_symbols": with_unfindable,
            "unfindable_builder_symbols": int(
                drift_df["builder_symbols_not_in_wikipathways"].sum()),
            "note": "Sizes compared on distinct Ensembl, not on labels: WikiPathways' "
                    "rdfs:label is not unique per gene and comparing label counts reports "
                    "drift on almost every pathway. A builder symbol absent from every "
                    "WikiPathways label for its pathway is the stale-corpus signal.",
        },
        "circularity_caveat": "Method 3's gene sets are the OUTPUT of the pathway mapping, not "
                              "evidence for it: pathways were selected partly on their genes. "
                              "Legitimate for the cross-method comparison, illegitimate as "
                              "corroboration of Method 3 itself (limitation 4).",
        "identifier_note": "No BridgeDb round-trip is involved. WikiPathways holds the Ensembl "
                           "cross-reference directly, so Method 3 has no symbol-to-Ensembl "
                           "mapping loss to report, unlike Methods 1 and 2.",
        "outputs": [
            "data/method3/derived/method3_pathway_genes.tsv",
            "data/method3/derived/method3_ke_ensembl_sets.tsv",
            "data/method3/derived/method3_genes_by_organ.tsv",
            "data/method3/derived/method3_corpus_drift.tsv",
        ],
    })


if __name__ == "__main__":
    main()
