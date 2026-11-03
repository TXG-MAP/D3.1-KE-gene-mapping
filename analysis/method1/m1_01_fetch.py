#!/usr/bin/env python3
"""Stage 1 — pull everything Method 1 needs out of the pinned AOP-Wiki RDF release.

Writes flat, long-format CSVs to ``data/method1/raw/``. Nothing downstream talks to the
endpoint, so every later stage is reproducible offline from these files.

The pin is verified before anything is fetched: the endpoint regenerates weekly and will
drift off commit 97191db. A drifted endpoint aborts the run unless ``--allow-drift`` is given.

Usage:
    python m1_01_fetch.py [--allow-drift]
"""

from __future__ import annotations

import argparse

import pandas as pd

import m1lib as m1


def _flag(series: pd.Series) -> pd.Series:
    """Coerce a SPARQL xsd:boolean to a real bool.

    Virtuoso serialises booleans as "1"/"0" rather than "true"/"false", and an unbound
    variable arrives as None. Accepting only "true" silently turns every flag into False —
    which is exactly the failure this helper exists to prevent, so it is asserted downstream.
    """
    return series.notna() & series.astype(str).str.lower().isin({"true", "1"})


def _hgnc_id(uri: pd.Series) -> pd.Series:
    return uri.str.rsplit("/", n=1).str[-1]


def fetch_ke_inventory() -> pd.DataFrame:
    m1.banner("KE inventory")
    df = m1.run_query(m1.read_query("method1c_ke_inventory.rq"))
    df = df.rename(columns={"ketitle": "title", "description_length": "description_length"})
    df["title"] = df["title"].str.strip()
    df["description_length"] = pd.to_numeric(df["description_length"], errors="coerce")
    df["has_description"] = df["description_length"].notna()
    # A KE with two declared biological levels would duplicate here; collapse deterministically
    # and record it rather than silently keeping both rows.
    dupes = df[df.duplicated("keid", keep=False)]
    if not dupes.empty:
        print(f"  note: {dupes['keid'].nunique()} KEs carry >1 biological level; "
              "keeping the first in the ordered vocabulary")
        order = {lvl: i for i, lvl in enumerate(m1.BIOLEVELS)}
        df["_o"] = df["biolevel"].map(order).fillna(99)
        df = df.sort_values(["keid", "_o"]).drop_duplicates("keid").drop(columns="_o")
    described = int(df["has_description"].sum())
    print(f"  {len(df)} Key Events, {df['biolevel'].notna().sum()} with a biological level")
    print(f"  {described} carry a free-text description ({described / len(df):.1%}) — this is "
          "the population Method 1 can actually read")
    return df.sort_values("keid").reset_index(drop=True)


def fetch_ke_genes() -> pd.DataFrame:
    m1.banner("KE-gene pairs with detector provenance")
    df = m1.run_query(m1.read_query("method1b_regex_vs_ner.rq"))
    df["by_regex"] = _flag(df["byRegex"])
    df["by_ner"] = _flag(df["byNer"])
    df = df.drop(columns=["byRegex", "byNer"])
    df["hgnc"] = _hgnc_id(df["gene"])
    df = df.rename(columns={"ketitle": "title"})

    # One row per (KE, gene). The per-gene OPTIONAL binding should already guarantee this;
    # assert it, because the cross-product failure mode is silent and would inflate everything.
    before = len(df)
    df = (df.groupby(["keid", "hgnc"], as_index=False)
            .agg(title=("title", "first"), gene=("gene", "first"),
                 symbol=("symbol", "first"),
                 by_regex=("by_regex", "any"), by_ner=("by_ner", "any")))
    if before != len(df):
        print(f"  collapsed {before} rows to {len(df)} unique (KE, gene) pairs")

    df["detector"] = pd.Series(
        pd.NA, index=df.index, dtype="object"
    ).mask(df.by_regex & df.by_ner, "both") \
     .mask(df.by_regex & ~df.by_ner, "regex only") \
     .mask(~df.by_regex & df.by_ner, "BERN2 only")

    orphan = int(df["detector"].isna().sum())
    if orphan == len(df):
        # Every pair unattributed means the flags failed to parse, not that the data changed.
        # This is what a boolean-serialisation change upstream looks like, and it would sail
        # through as "BERN2 contributed nothing" if it were merely warned about.
        raise RuntimeError(
            "no (KE, gene) pair carries a detector flag — the provenance predicates almost "
            "certainly failed to parse rather than being absent from the data. Check how the "
            "endpoint serialises xsd:boolean (Virtuoso returns 1/0, not true/false) before "
            "trusting any regex-vs-BERN2 number."
        )
    if orphan:
        # A gene in edam:data_1025 claimed by neither provenance predicate: the union and the
        # provenance layer disagree. Real if it happens, so it is surfaced, not swallowed.
        print(f"  WARNING: {orphan} (KE, gene) pairs carry neither provenance predicate")
        df["detector"] = df["detector"].fillna("unattributed")

    print(f"  {len(df)} KE-gene pairs · {df['keid'].nunique()} KEs · "
          f"{df['hgnc'].nunique()} unique genes")
    print("  " + " · ".join(f"{k}: {v}" for k, v in df["detector"].value_counts().items()))
    return df.sort_values(["keid", "hgnc"]).reset_index(drop=True)


def fetch_aop_ke() -> tuple[pd.DataFrame, int]:
    m1.banner("AOP -> KE membership")
    total_aops = int(m1.run_query(
        "PREFIX aopo: <http://aopkb.org/aop_ontology#> "
        "SELECT (COUNT(DISTINCT ?a) AS ?n) WHERE { ?a a aopo:AdverseOutcomePathway }"
    )["n"].iloc[0])

    df = m1.run_query(m1.read_query("network_organ_aops.rq"))
    df = df.rename(columns={"aoptitle": "aop_title"})
    df["aop_number"] = df["aopid"].str.extract(r"(\d+)").astype(int)
    with_kes = int(df["aopid"].nunique())
    print(f"  {total_aops} AOPs declared · {with_kes} of them carry >=1 Key Event · "
          f"{total_aops - with_kes} are empty")
    print(f"  {df['keid'].nunique()} distinct KEs · {len(df)} (AOP, KE) memberships")
    print("  roles: " + " · ".join(f"{k}: {v}" for k, v in df["role"].value_counts().items()))
    return df.sort_values(["aop_number", "keid"]).reset_index(drop=True), total_aops


def fetch_ker_genes() -> pd.DataFrame:
    m1.banner("KER-gene pairs (supplementary)")
    df = m1.run_query(m1.read_query("method1d_ker_gene_annotations.rq"))
    df["by_regex"] = _flag(df["byRegex"])
    df["by_ner"] = _flag(df["byNer"])
    df = df.drop(columns=["byRegex", "byNer"])
    df["hgnc"] = _hgnc_id(df["gene"])
    df = (df.groupby(["kerid", "hgnc"], as_index=False)
            .agg(upstream=("upstream", "first"), downstream=("downstream", "first"),
                 gene=("gene", "first"),
                 by_regex=("by_regex", "any"), by_ner=("by_ner", "any")))
    print(f"  {df['kerid'].nunique()} KERs carry genes · {len(df)} KER-gene pairs · "
          f"{df['hgnc'].nunique()} unique genes")
    return df.sort_values(["kerid", "hgnc"]).reset_index(drop=True)


def fetch_gene_xrefs(all_genes: set[str]) -> pd.DataFrame:
    m1.banner("gene cross-references (route A: owl:sameAs in the pinned release)")
    df = m1.run_query(m1.read_query("method1e_gene_xrefs.rq"))
    df["hgnc"] = _hgnc_id(df["gene"])
    df["space"] = df["xref"].str.extract(r"identifiers\.org/([^/]+)/")
    df["xref_id"] = df["xref"].str.rsplit("/", n=1).str[-1]
    wide = df[["hgnc", "symbol", "space", "xref_id"]].drop_duplicates()
    # The denominator is every gene Method 1 annotates onto a KE, NOT the genes this query
    # returned — a gene with no cross-reference at all simply does not appear here, and using
    # the returned set as the denominator would report a 100 % mapping rate by construction.
    for space in ("ensembl", "ncbigene"):
        mapped = set(wide.loc[wide.space == space, "hgnc"])
        print(f"  {space}: {len(mapped & all_genes)} of {len(all_genes)} annotated genes "
              f"have >=1 cross-reference ({len(all_genes - mapped)} unmapped)")
    return wide.sort_values(["hgnc", "space", "xref_id"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-drift", action="store_true",
                        help="proceed even if the endpoint no longer serves the pinned release")
    args = parser.parse_args()

    m1.banner("verifying the pin")
    pin_status = m1.verify_pin(allow_drift=args.allow_drift)

    ke = fetch_ke_inventory()
    ke_genes = fetch_ke_genes()
    aop_ke, total_aops = fetch_aop_ke()
    ker_genes = fetch_ker_genes()
    xrefs = fetch_gene_xrefs(set(ke_genes["hgnc"]))

    outputs = {
        "ke_inventory.csv": ke,
        "ke_gene_long.csv": ke_genes,
        "aop_ke.csv": aop_ke,
        "ker_gene_long.csv": ker_genes,
        "gene_xrefs.csv": xrefs,
    }
    m1.banner("writing raw tables")
    written = {}
    for name, frame in outputs.items():
        path = m1.save_table(frame, m1.RAW_DIR / name)
        written[name] = {"rows": len(frame), "sha256": m1.sha256(path)}
        print(f"  -> {path.relative_to(m1.ROOT)}  ({len(frame)} rows)")

    # Cross-check the fetched counts against the manifest declared in pinned-release.yaml.
    # This is the sanity gate the BridgeDb incident of 2026-07-25 argues for: a silent
    # upstream format change produced an HTTP 200 and a 78.7 % drop in gene counts.
    m1.banner("sanity check against pinned-release.yaml manifest")
    declared = m1.PIN["manifest"]
    observed = {
        "aops": total_aops,
        "key_events": int(len(ke)),
        "kes_with_genes": int(ke_genes["keid"].nunique()),
        "kes_gene_by_ner": int(ke_genes.loc[ke_genes.by_ner, "keid"].nunique()),
        "kes_gene_by_regex": int(ke_genes.loc[ke_genes.by_regex, "keid"].nunique()),
        "unique_genes_on_kes": int(ke_genes["hgnc"].nunique()),
    }
    # Method 1 mines the KE description, so a KE without one cannot be annotated. Assert the
    # containment rather than assuming it: if a gene-annotated KE ever appears without a
    # description, the upstream method has started reading some other field and every
    # "coverage among mineable KEs" figure in the report silently changes meaning.
    described = set(ke.loc[ke.has_description, "keid"])
    annotated_without_text = sorted(set(ke_genes["keid"]) - described)
    if annotated_without_text:
        print(f"\n  WARNING: {len(annotated_without_text)} gene-annotated KEs carry no "
              f"dc:description: {annotated_without_text[:10]}")
        print("  Method 1 is no longer confined to the description field. The 'coverage among "
              "KEs with a description' denominator must be re-derived before it is quoted.")

    extra = {
        "kes_with_description": len(described),
        "gene_annotated_kes_without_description": len(annotated_without_text),
        "aops_with_at_least_one_ke": int(aop_ke["aopid"].nunique()),
        "aops_with_no_key_event": total_aops - int(aop_ke["aopid"].nunique()),
        "kers_with_genes": int(ker_genes["kerid"].nunique()),
        "unique_genes_on_kers": int(ker_genes["hgnc"].nunique()),
    }
    mismatches = []
    for key, value in observed.items():
        exp = declared.get(key)
        flag = "ok" if exp == value else f"MISMATCH (manifest says {exp})"
        if exp != value:
            mismatches.append((key, exp, value))
        print(f"  {key:24s} {value:>6d}  {flag}")
    for key, value in extra.items():
        print(f"  {key:24s} {value:>6d}  (not in the manifest)")
    if mismatches:
        print("\n  The manifest in pinned-release.yaml disagrees with the live endpoint on "
              f"{len(mismatches)} count(s). Small differences are usually a counting "
              "definition, not data loss — resolve before quoting either number.")

    m1.write_manifest("01_fetch", {
        "pin_verification": pin_status,
        "outputs": written,
        "observed_counts": {**observed, **extra},
        "declared_counts": declared,
        "count_mismatches": [
            {"key": k, "declared": e, "observed": o} for k, e, o in mismatches
        ],
        "queries": ["method1b_regex_vs_ner.rq", "method1c_ke_inventory.rq",
                    "network_organ_aops.rq", "method1d_ker_gene_annotations.rq",
                    "method1e_gene_xrefs.rq"],
    })
    print(f"\nStage 1 complete. Raw tables in {m1.RAW_DIR.relative_to(m1.ROOT)}/")


if __name__ == "__main__":
    main()
