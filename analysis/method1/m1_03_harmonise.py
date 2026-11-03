#!/usr/bin/env python3
"""Stage 3 — harmonise Method 1's HGNC gene identifiers to Ensembl (and NCBI Gene).

The three-method comparison in Section 07 is run in Ensembl gene IDs,
because Method 2 already produces them and transcriptomics data arrives in them. Method 1
natively carries HGNC, so it needs a conversion — and the conversion loses genes, which is a
result the deliverable is asked to report rather than a bookkeeping detail.

Two routes to Ensembl exist and they do not agree:

``rdf``       the ``owl:sameAs`` cross-references already in the pinned release. Frozen with
              the pin, works offline, and is what the report quotes — it is the only route the
              pin covers.
``bridgedb``  a live call to the BridgeDb webservice. Current, but unpinned and therefore not
              reproducible from the pin alone. Run as a cross-check.

The disagreement between them is the interesting number: it bounds how much of Method 1's
gene set is sensitive to *when* the identifier mapping was performed.

**Why this stage checks itself.** On 2026-07-25 BridgeDb 2.1.8 changed its batch response from
TSV to JSON and the AOP-Wiki RDF pipeline silently produced a 78.7 % drop in gene counts on an
HTTP 200. So the BridgeDb route asserts a floor on its own mapping rate and refuses to write a
result that looks like that failure.

Usage:
    python m1_03_harmonise.py [--skip-bridgedb] [--chunk 500]
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import defaultdict

import pandas as pd
import requests

import m1lib as m1

BRIDGEDB = "https://webservice.bridgedb.org"
SPECIES = "Human"

#: Below this fraction of genes mapping to Ensembl, assume the service or its response format
#: broke rather than that the data changed. The RDF route resolves ~99.7 %, and the July 2026
#: incident showed up as a drop to roughly a fifth of the expected count.
MAPPING_RATE_FLOOR = 0.80


# --------------------------------------------------------------------------------------
# Route A — the pinned RDF's own cross-references
# --------------------------------------------------------------------------------------


def harmonise_via_rdf(genes: set[str], xrefs: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    m1.banner("route A — owl:sameAs cross-references in the pinned release")
    rows = []
    stats = {}
    for space in ("ensembl", "ncbigene"):
        sub = xrefs[(xrefs.space == space) & (xrefs.hgnc.isin(genes))]
        mapped = set(sub["hgnc"])
        rows.append(sub[["hgnc", "space", "xref_id"]])
        multi = sub.groupby("hgnc").size()
        stats[space] = {
            "genes_in": len(genes),
            "genes_mapped": len(mapped),
            "genes_unmapped": len(genes - mapped),
            "mapping_rate": len(mapped) / len(genes),
            "target_ids": int(sub["xref_id"].nunique()),
            "genes_with_multiple_targets": int((multi > 1).sum()),
            "max_targets_for_one_gene": int(multi.max()) if len(multi) else 0,
        }
        print(f"  {space:9s} {len(mapped):5d}/{len(genes)} genes mapped "
              f"({len(mapped) / len(genes):.1%}) -> {sub['xref_id'].nunique()} distinct IDs · "
              f"{int((multi > 1).sum())} genes map to >1 target")
    return pd.concat(rows, ignore_index=True), stats


# --------------------------------------------------------------------------------------
# Route B — the BridgeDb webservice
# --------------------------------------------------------------------------------------


def bridgedb_version() -> dict:
    """Capture what the service is serving, so the run can be described even though the
    service itself cannot be pinned."""
    try:
        text = requests.get(f"{BRIDGEDB}/{SPECIES}/properties", timeout=30).text
    except Exception as exc:
        print(f"  ! could not read BridgeDb properties: {exc}", file=sys.stderr)
        return {"error": str(exc)}
    props, current = [], {}
    for line in text.splitlines():
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        if key == "DATASOURCENAME" and current:
            props.append(current)
            current = {}
        current[key] = value
    if current:
        props.append(current)
    ensembl = next((p for p in props if p.get("DATASOURCENAME") == "Ensembl"), {})
    info = {
        "ensembl_datasource_version": ensembl.get("DATASOURCEVERSION"),
        "ensembl_build_date": ensembl.get("BUILDDATE"),
        "queried_on": pd.Timestamp.now("UTC").isoformat(),
    }
    print(f"  BridgeDb Ensembl datasource v{info['ensembl_datasource_version']} "
          f"(built {info['ensembl_build_date']})")
    return info


def bridgedb_batch(hgnc_ids: list[str], target: str, chunk: int) -> dict[str, set[str]]:
    """Map HGNC accession numbers to one target system through the xrefsBatch endpoint.

    Input identifiers use the ``Hac`` system code (HGNC Accession number, ``HGNC:1234``); the
    plain ``H`` code is the gene *symbol*, which is ambiguous and drifts between releases.

    Two things about this endpoint's response are easy to get wrong, and getting either wrong
    inflates the mapping rate to ~100 % instead of erroring:

    1. ``?dataSource=`` is honoured by the single-identifier ``/xrefs/`` endpoint but **ignored
       by the batch endpoint**, which always returns every system it knows.
    2. Batch results are therefore prefixed with the system code (``En:ENSG…``, ``L:5468``,
       ``Pd:``, ``T:GO:…``), whereas single lookups return bare identifiers.

    So the target system is selected here, by prefix, rather than being requested from the
    service. Note ``L:`` must be matched exactly — a naive ``startswith("L")`` also swallows
    ``Il:`` (Illumina) and would silently mix probe IDs into the NCBI Gene set.
    """
    out: dict[str, set[str]] = defaultdict(set)
    url = f"{BRIDGEDB}/{SPECIES}/xrefsBatch/Hac"
    prefix = f"{target}:"
    total = len(hgnc_ids)
    for start in range(0, total, chunk):
        block = hgnc_ids[start:start + chunk]
        payload = "\n".join(f"HGNC:{h}" for h in block)
        for attempt in range(1, 4):
            try:
                resp = requests.post(url, data=payload.encode(), timeout=180)
                resp.raise_for_status()
                break
            except Exception as exc:
                if attempt == 3:
                    raise RuntimeError(f"BridgeDb batch failed: {exc}") from exc
                print(f"    ! attempt {attempt} failed ({exc}); retrying", file=sys.stderr)
                time.sleep(3 * attempt)
        for line in resp.text.splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            key, mapped = parts[0], parts[2]
            hgnc = key.split(":", 1)[-1]
            if not mapped or mapped == "N/A":
                continue
            hits = {m[len(prefix):] for m in mapped.split(",") if m.startswith(prefix)}
            if hits:
                out[hgnc].update(hits)
        print(f"    {min(start + chunk, total):5d}/{total} submitted", end="\r")
    print(" " * 48, end="\r")
    return out


def harmonise_via_bridgedb(genes: set[str], chunk: int) -> tuple[pd.DataFrame, dict]:
    m1.banner("route B — live BridgeDb webservice")
    version = bridgedb_version()
    ordered = sorted(genes, key=int)
    rows, stats = [], {"service": version}
    for space, code in (("ensembl", "En"), ("ncbigene", "L")):
        mapping = bridgedb_batch(ordered, code, chunk)
        rate = len(mapping) / len(genes)
        stats[space] = {
            "genes_in": len(genes),
            "genes_mapped": len(mapping),
            "genes_unmapped": len(genes) - len(mapping),
            "mapping_rate": rate,
            "target_ids": len({v for vals in mapping.values() for v in vals}),
        }
        print(f"  {space:9s} {len(mapping):5d}/{len(genes)} genes mapped ({rate:.1%})")
        if rate < MAPPING_RATE_FLOOR:
            raise RuntimeError(
                f"BridgeDb mapped only {rate:.1%} of genes to {space}, below the "
                f"{MAPPING_RATE_FLOOR:.0%} floor. This is the signature of a response-format "
                "change rather than of the data — the July 2026 incident looked exactly like "
                "this and returned HTTP 200 throughout. Inspect a raw response before "
                "trusting or lowering this floor."
            )
        rows.append(pd.DataFrame(
            [{"hgnc": h, "space": space, "xref_id": t} for h, ts in mapping.items() for t in ts]
        ))
    return pd.concat(rows, ignore_index=True), stats


# --------------------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------------------


def compare_routes(rdf: pd.DataFrame, bridge: pd.DataFrame, genes: set[str]) -> pd.DataFrame:
    m1.banner("route agreement")
    records = []
    for space in ("ensembl", "ncbigene"):
        a = set(rdf.loc[rdf.space == space, "hgnc"])
        b = set(bridge.loc[bridge.space == space, "hgnc"])
        pa = set(map(tuple, rdf.loc[rdf.space == space, ["hgnc", "xref_id"]].values))
        pb = set(map(tuple, bridge.loc[bridge.space == space, ["hgnc", "xref_id"]].values))
        records.append({
            "target": space,
            "genes mapped by RDF": len(a),
            "genes mapped by BridgeDb": len(b),
            "mapped by both": len(a & b),
            "RDF only": len(a - b),
            "BridgeDb only": len(b - a),
            "mapped by neither": len(genes - a - b),
            "identical pairs": len(pa & pb),
            "pair agreement": len(pa & pb) / len(pa | pb) if (pa | pb) else float("nan"),
        })
    out = pd.DataFrame(records)
    print(m1.markdown_table(out, floatfmt="{:.3f}"))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-bridgedb", action="store_true",
                        help="use only the pinned RDF cross-references (fully offline)")
    parser.add_argument("--chunk", type=int, default=500,
                        help="genes per BridgeDb batch request (default 500)")
    args = parser.parse_args()

    ke_genes = pd.read_csv(m1.RAW_DIR / "ke_gene_long.csv", dtype={"hgnc": str})
    xrefs = pd.read_csv(m1.RAW_DIR / "gene_xrefs.csv", dtype={"hgnc": str, "xref_id": str})
    genes = set(ke_genes["hgnc"])
    print(f"Harmonising {len(genes)} unique HGNC genes annotated onto "
          f"{ke_genes['keid'].nunique()} Key Events.")

    rdf_map, rdf_stats = harmonise_via_rdf(genes, xrefs)
    rdf_map["route"] = "rdf"

    parts, bridge_stats, comparison = [rdf_map], None, pd.DataFrame()
    if not args.skip_bridgedb:
        try:
            bridge_map, bridge_stats = harmonise_via_bridgedb(genes, args.chunk)
            bridge_map["route"] = "bridgedb"
            parts.append(bridge_map)
            comparison = compare_routes(rdf_map, bridge_map, genes)
        except Exception as exc:
            print(f"\n  BridgeDb cross-check could not be completed: {exc}", file=sys.stderr)
            print("  Continuing with the pinned RDF route only, which is what the report "
                  "quotes.", file=sys.stderr)
            bridge_stats = {"error": str(exc)}

    combined = pd.concat(parts, ignore_index=True).drop_duplicates()
    m1.save_table(combined, m1.DERIVED_DIR / "gene_id_harmonisation.csv")
    print(f"\n  -> {(m1.DERIVED_DIR / 'gene_id_harmonisation.csv').relative_to(m1.ROOT)}")

    # The per-KE Ensembl gene sets the comparison section will consume, on the route the
    # report quotes. Genes with no Ensembl equivalent are dropped here — that is the loss.
    ens = combined[(combined.route == "rdf") & (combined.space == "ensembl")]
    per_ke = (ke_genes[["keid", "hgnc"]]
              .merge(ens[["hgnc", "xref_id"]], on="hgnc", how="left")
              .rename(columns={"xref_id": "ensembl"}))
    lost_pairs = int(per_ke["ensembl"].isna().sum())
    kes_losing_all = sorted(set(per_ke.loc[per_ke.ensembl.isna(), "keid"])
                            - set(per_ke.loc[per_ke.ensembl.notna(), "keid"]))
    ensembl_sets = (per_ke.dropna(subset=["ensembl"])
                    .groupby("keid")["ensembl"].apply(lambda s: sorted(set(s)))
                    .reset_index())
    ensembl_sets["n_ensembl"] = ensembl_sets["ensembl"].str.len()
    ensembl_sets["ensembl"] = ensembl_sets["ensembl"].str.join(", ")
    m1.save_table(ensembl_sets, m1.DERIVED_DIR / "method1_ke_ensembl_sets.csv")

    print(f"  {lost_pairs} KE-gene pairs have no Ensembl equivalent; "
          f"{len(kes_losing_all)} KEs lose their entire gene set")
    if kes_losing_all:
        print(f"    {kes_losing_all}")
    print(f"  -> {(m1.DERIVED_DIR / 'method1_ke_ensembl_sets.csv').relative_to(m1.ROOT)} "
          f"({len(ensembl_sets)} KEs)")

    if not comparison.empty:
        m1.save_table(comparison, m1.DERIVED_DIR / "harmonisation_route_agreement.csv")

    m1.write_manifest("03_harmonise", {
        "genes_in": len(genes),
        "route_rdf": rdf_stats,
        "route_bridgedb": bridge_stats,
        "quoted_route": "rdf",
        "quoted_route_rationale": "the only route the release pin covers; BridgeDb is live "
                                  "and therefore not reproducible from the pin alone",
        "ke_gene_pairs_without_ensembl": lost_pairs,
        "kes_losing_entire_gene_set": kes_losing_all,
        "mapping_rate_floor": MAPPING_RATE_FLOOR,
    })
    print("\nStage 3 complete.")


if __name__ == "__main__":
    main()
