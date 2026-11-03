#!/usr/bin/env python3
"""Methods 2 and 3 gene sets, taken from the molAOP Builder itself.

**This is the authoritative route for both methods** (decided by Marvin, 2026-08-14). Methods 2
and 3 were both curated *in* the builder, and the deliverable reports what the resource actually
holds — not what an independent re-derivation says it should hold. Taking both from the same
service also makes them symmetric: same provenance, same corpus snapshot, same export machinery,
and one identifier round-trip instead of two different ones.

    /exports/gmt/ke-go   KE -> GO Biological Process, 336 approved pairs   (Method 2)
    /exports/gmt/ke-wp   KE -> WikiPathways,          459 approved pairs   (Method 3)

Both exports carry HGNC **symbols**, so both need the same symbol -> Ensembl step, and both
therefore have a mapping loss to report under Q2. That is a change from the earlier
independent derivations, and a fair one: the loss is real, it is incurred by anyone consuming
these exports, and hiding it by going around the resource would have understated it.

**What this supersedes, and what it does not.** The independent derivations still run and are
still worth keeping:

  * ``analysis/method3/m3_13_genesets.py`` reads pathway membership from WikiPathways SPARQL.
  * ``analysis/method2/m2_02_genes.py`` expands GO terms through a pinned GOA with its own
    closure.

They are now **cross-checks**, not the reported numbers. Their value is that they measure the
drift between what the builder serves and what the upstream resources currently say, which is
exactly the question upstream molAOP-builder #263 is open about. Keeping both also means a
future corpus refresh can be quantified rather than discovered.

The builder pins its own upstream snapshot and exposes it at ``/stats#source-versions``. Those
versions are read and recorded here, because a gene set from this route is only interpretable
against them — and they are NOT the same as the releases the independent routes use.

Network: 1 GET for /stats, 2 GETs for the exports, then BridgeDb in batches.
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from collections import defaultdict
from datetime import date
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent / "method3"))
import m3lib as L  # noqa: E402

ROOT = L.ROOT
STAMP = date.today().isoformat()

BRIDGEDB = "https://webservice.bridgedb.org"
SPECIES = "Homo sapiens"
MAPPING_RATE_FLOOR = 0.80

#: resource -> (export path, id pattern, method, output directory)
RESOURCES = {
    "ke-go": {
        "path": "/exports/gmt/ke-go",
        "pattern": re.compile(r"^KE(\d+)_.*_(GO:\d+)$"),
        "method": "method2",
        "out": ROOT / "data" / "method2" / "derived",
        "corpus": ROOT / "data" / "method2" / "corpus",
        "target": "GO Biological Process term",
    },
    "ke-wp": {
        "path": "/exports/gmt/ke-wp",
        "pattern": re.compile(r"^KE(\d+)_.*_(WP\d+)$"),
        "method": "method3",
        "out": ROOT / "data" / "method3" / "derived",
        "corpus": ROOT / "data" / "method3" / "corpus",
        "target": "WikiPathways pathway",
    },
}

FOUR = ["brain", "kidney", "liver", "lung"]
NODE_FILE = {"brain": "brain_nodes.tsv", "kidney": "kidney_nodes.tsv",
             "liver": "liver_nodes.tsv", "lung": "lungs_nodes.tsv"}


def builder_source_versions() -> dict:
    """The upstream snapshot the running container serves, from its own footer.

    Without these a gene set from this route cannot be interpreted: the builder's WikiPathways
    snapshot is a month behind the live endpoint, and its GO snapshot two months behind the GOA
    the independent route uses. That is not a defect to correct here — it is the pin, and it has
    to be reported as one.
    """
    html = requests.get(f"{L.BASE}/stats", timeout=90,
                        headers={"User-Agent": L.USER_AGENT}).text
    out = {}
    for label, key in (("WP", "wikipathways"), ("GO", "go"),
                       ("Reactome", "reactome"), ("AOP-Wiki", "aopwiki")):
        m = re.search(rf"{re.escape(label)}&nbsp;([^\s<&]+)", html)
        if m:
            out[key] = m.group(1)
    if not out:
        raise SystemExit(
            "could not read the source-version strip from /stats. Refusing to continue: a gene "
            "set taken from the builder is only interpretable against the corpus snapshot it "
            "was built from."
        )
    return out


def fetch_export(resource: str) -> tuple[list[str], str]:
    spec = RESOURCES[resource]
    body, headers = L.get(spec["path"], raw=True, pace=1.0)
    # Each export is snapshotted under ITS OWN method, not under whichever method's library
    # this script happens to import for its HTTP plumbing.
    spec["corpus"].mkdir(parents=True, exist_ok=True)
    dest = spec["corpus"] / f"{resource}_builder_{STAMP}.gmt"
    dest.write_bytes(body)
    lines = [ln for ln in body.decode().splitlines()
             if ln.strip() and not ln.startswith("#")]
    return lines, headers.get("x-export-revision", "")


def parse_gmt(lines: list[str], pattern: re.Pattern) -> pd.DataFrame:
    """One row per KE-target pair, with the symbols the builder serves for it.

    Every field is stripped: the last gene on a line keeps its terminator, and an unstripped
    symbol matches nothing downstream while looking like a real mapping loss.
    """
    rows, unparsed = [], []
    for line in lines:
        parts = [p.strip() for p in line.split("\t")]
        m = pattern.match(parts[0])
        if not m:
            unparsed.append(parts[0])
            continue
        rows.append({
            "ke_id": f"KE {m.group(1)}",
            "target_id": m.group(2),
            "target_title": parts[1] if len(parts) > 1 else "",
            "symbols": sorted({g for g in parts[2:] if g}),
        })
    if unparsed:
        raise SystemExit(
            f"{len(unparsed)} GMT set names did not parse, e.g. {unparsed[:3]}. The name format "
            f"is the only join key this export offers, so a silent skip would drop mappings."
        )
    return pd.DataFrame(rows)


def bridgedb_symbols(symbols: list[str], chunk: int = 500) -> dict[str, set[str]]:
    """HGNC symbol -> Ensembl. The batch endpoint IGNORES ?dataSource and returns every system
    it knows, prefixed, so the target is selected here by prefix rather than requested."""
    out: dict[str, set[str]] = defaultdict(set)
    url = f"{BRIDGEDB}/{SPECIES}/xrefsBatch/H"
    for start in range(0, len(symbols), chunk):
        block = symbols[start:start + chunk]
        for attempt in range(1, 4):
            try:
                resp = requests.post(url, data="\n".join(block).encode(), timeout=180)
                resp.raise_for_status()
                break
            except Exception as exc:
                if attempt == 3:
                    raise RuntimeError(f"BridgeDb batch failed: {exc}") from exc
                time.sleep(3 * attempt)
        for line in resp.text.splitlines():
            parts = line.split("\t")
            if len(parts) < 3 or not parts[2] or parts[2] == "N/A":
                continue
            hits = {m[3:] for m in parts[2].split(",") if m.startswith("En:")}
            if hits:
                out[parts[0]].update(hits)
        print(f"    {min(start + chunk, len(symbols)):6d}/{len(symbols)} submitted", end="\r")
    print(" " * 48, end="\r")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-bridgedb", action="store_true")
    args = ap.parse_args()

    L.banner("Methods 2 and 3 gene sets, from the molAOP Builder")

    versions = builder_source_versions()
    print("  builder corpus snapshot: "
          + " · ".join(f"{k} {v}" for k, v in versions.items()))

    networks = {o: set(pd.read_csv(ROOT / "data" / "networks" / NODE_FILE[o],
                                   sep="\t")["id"].str.strip()) for o in FOUR}
    union_kes = set().union(*networks.values())

    pairs, revisions = {}, {}
    for resource, spec in RESOURCES.items():
        lines, revision = fetch_export(resource)
        df = parse_gmt(lines, spec["pattern"])
        pairs[resource] = df
        revisions[resource] = revision
        print(f"  {resource:<7} {len(df):4d} approved pairs over {df['ke_id'].nunique():4d} KEs, "
              f"{df['target_id'].nunique():4d} distinct {spec['target']}s "
              f"(revision {revision})")

    # --- one symbol pool, one round-trip ------------------------------------------------
    all_symbols = sorted({s for df in pairs.values() for row in df["symbols"] for s in row})
    print(f"\n  {len(all_symbols)} distinct symbols across both exports")

    ens: dict[str, set[str]] = {}
    bridge_stats = {}
    if not args.skip_bridgedb:
        print("  harmonising symbols to Ensembl via BridgeDb")
        ens = bridgedb_symbols(all_symbols)
        rate = len(ens) / len(all_symbols)
        print(f"    {len(ens)}/{len(all_symbols)} symbols mapped ({rate:.1%}), "
              f"{len({e for v in ens.values() for e in v})} distinct Ensembl genes")
        if rate < MAPPING_RATE_FLOOR:
            raise SystemExit(
                f"BridgeDb mapped only {rate:.1%} of symbols to Ensembl, below the "
                f"{MAPPING_RATE_FLOOR:.0%} floor — the signature of a response-format change "
                f"rather than of the data."
            )
        unmapped = sorted(set(all_symbols) - set(ens))
        bridge_stats = {
            "symbols_in": len(all_symbols),
            "symbols_mapped": len(ens),
            "mapping_rate": rate,
            "ensembl_genes": len({e for v in ens.values() for e in v}),
            "unmapped_examples": unmapped[:25],
        }
        pd.DataFrame({"symbol": unmapped}).to_csv(
            ROOT / "data" / "comparison" / "builder_symbols_unmapped.tsv",
            sep="\t", index=False)
        print(f"    {len(unmapped)} symbols carry no Ensembl cross-reference "
              f"(written to data/comparison/builder_symbols_unmapped.tsv)")

    # --- per-method KE gene sets --------------------------------------------------------
    summary = {}
    for resource, spec in RESOURCES.items():
        df = pairs[resource]
        rows = []
        for ke, grp in df.groupby("ke_id"):
            symbols = sorted({s for row in grp["symbols"] for s in row})
            genes = sorted({e for s in symbols for e in ens.get(s, ())}) if ens else []
            rows.append({
                "ke_id": ke,
                "n_targets": grp["target_id"].nunique(),
                "targets": "|".join(sorted(grp["target_id"])),
                "n_symbols": len(symbols),
                "n_genes": len(genes),
                "symbols_without_ensembl": len([s for s in symbols if s not in ens]),
                "ensembl": "|".join(genes),
            })
        ke_sets = pd.DataFrame(rows).sort_values(
            "ke_id", key=lambda s: s.str.extract(r"(\d+)")[0].astype(int))
        out = spec["out"] / f"{spec['method']}_builder_ke_ensembl_sets.tsv"
        L.save_table(ke_sets, out)

        in_frame = ke_sets[ke_sets["ke_id"].isin(union_kes)]
        organ_rows = []
        for organ in FOUR + ["four-organ union"]:
            ids = union_kes if organ == "four-organ union" else networks[organ]
            sub = ke_sets[ke_sets["ke_id"].isin(ids)]
            genes = set("|".join(sub["ensembl"]).split("|")) - {""} if len(sub) else set()
            with_genes = int((sub["n_genes"] > 0).sum()) if len(sub) else 0
            organ_rows.append({
                "organ": organ,
                "kes_in_frame": len(ids),
                "kes_with_a_mapping": len(sub),
                "kes_with_genes": with_genes,
                "pct_of_frame_with_genes": round(100 * with_genes / len(ids), 1),
                "unique_genes": len(genes),
                "ke_gene_pairs": int(sub["n_genes"].sum()) if len(sub) else 0,
                "median_genes_per_annotated_ke": float(
                    sub.loc[sub["n_genes"] > 0, "n_genes"].median()) if with_genes else 0.0,
            })
        organs = pd.DataFrame(organ_rows)
        L.save_table(organs, spec["out"] / f"{spec['method']}_builder_genes_by_organ.tsv")

        u = organs[organs["organ"] == "four-organ union"].iloc[0]
        print(f"    {spec['method']}: {len(ke_sets)} KEs service-wide, "
              f"{len(in_frame)} in the four networks, "
              f"{u['kes_with_genes']} with genes ({u['pct_of_frame_with_genes']} % of 538), "
              f"{u['unique_genes']} genes")
        summary[spec["method"]] = {
            "resource": resource,
            "export_revision": revisions[resource],
            "approved_pairs": len(df),
            "kes_service_wide": int(df["ke_id"].nunique()),
            "kes_in_the_four_networks": len(in_frame),
            "by_organ": organs.to_dict("records"),
        }

    write_manifest({
        "purpose": "Methods 2 and 3 gene sets taken from the molAOP Builder's own exports",
        "decision": "Marvin, 2026-08-14 — both methods were curated in the builder, so the "
                    "deliverable reports what the resource holds rather than an independent "
                    "re-derivation. This also makes the two methods symmetric: same corpus "
                    "snapshot, same export machinery, same identifier round-trip.",
        "builder_source_versions": versions,
        "builder_source_versions_note": "read from /stats#source-versions. These are NOT the "
                                        "releases the independent cross-checks use — the "
                                        "builder's WikiPathways snapshot is roughly a month "
                                        "behind the live endpoint and its GO snapshot about "
                                        "two months behind the GOA file. That gap is the "
                                        "subject of upstream molAOP-builder #263 and must be "
                                        "reported, not corrected here.",
        "exports": {r: revisions[r] for r in RESOURCES},
        "bridgedb": bridge_stats,
        "identifier_note": "both exports carry HGNC symbols, so both methods now incur the same "
                           "symbol-to-Ensembl loss. Earlier derivations avoided it — Method 3 by "
                           "reading Ensembl straight from WikiPathways, Method 2 by joining on "
                           "UniProt — but that understated what a consumer of these exports "
                           "actually gets.",
        "supersedes": [
            "data/method3/derived/method3_ke_ensembl_sets.tsv (WikiPathways SPARQL route)",
            "data/method2/derived/method2_ke_gene_sets.tsv (pinned GOA route)",
        ],
        "supersedes_note": "both remain in the repository and still run. They are cross-checks "
                           "that measure builder-versus-upstream drift, not the reported "
                           "numbers.",
        "methods": summary,
    })


def write_manifest(payload: dict) -> None:
    """Written to data/comparison/, because this stage spans two methods and belongs to
    neither."""
    import json
    from datetime import datetime, timezone
    payload = dict(payload)
    payload.update({
        "stage": "builder_genesets",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "aopwiki_rdf_commit": L.COMMIT,
        "molaop_builder_base": L.BASE,
        "d31_repo_revision": L.git_revision(),
        "python": sys.version.split()[0],
    })
    path = ROOT / "data" / "comparison" / "manifest_builder_genesets.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print(f"  -> {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
