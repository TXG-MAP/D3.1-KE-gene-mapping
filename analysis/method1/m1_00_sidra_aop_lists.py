#!/usr/bin/env python3
"""One-off — write out the four organ AOP lists that define Sidra's networks.

These four lists are the frozen definition of the organ AOP networks for D3.1 (decision,
2026-08-06). Method 2 was run on them and Method 1 is re-run onto them, so the two methods
sit on the same networks.

Where each list comes from:

``liver``
    Directly from ``data/networks/kes.tsv``, which carries an ``aopNumber`` column. 36 AOPs.

``brain``, ``kidney``, ``lungs``
    Her committed node tables carry no AOP column, so the AOP list is recovered: run her
    keyword query for the organ against the pinned release to get the candidate AOPs, then
    keep every candidate all of whose Key Events appear in her node table. The recovery is
    checked, not assumed — the union of the kept AOPs' KEs must equal her node table exactly,
    and the script aborts if it does not. All three reconcile exactly (203 / 68 / 134 KEs).

Liver is deliberately not recovered this way. Two of its KEs (KE 1419, KE 1420 in AOP 232)
were replaced in AOP-Wiki after her run, so the check would fail on drift rather than on a
selection difference.

Usage:
    python m1_00_sidra_aop_lists.py
"""

from __future__ import annotations

import csv
import json
import urllib.parse
import urllib.request

import m1lib as m1

NETWORKS = m1.ROOT / "data" / "networks"

# Verbatim from data/networks/SPARQL_queries.txt, except liver: the LIVER block in that file
# carries the brain keyword list by mistake, so liver is taken from AOP_selection.R.
KEYWORDS = {
    "brain": ("brain|neuro|central nervous system|cns|neurodegeneration|neuropathology|"
              "neurogenesis|neurotransmitter|neurodevelopment|neuroinflammation|"
              "blood brain barrier|bbb|cognitive|memory|learning|synaptic|glial|glia"),
    "kidney": "kidney|renal|nephro|glomerulus|glomerular|urinary|uremia",
    "lungs": ("lung|pulmonary|respiratory|bronchi|alveol|hypoxia|dyspnea|emphysema|"
              "pulmonary surfactant|asthma|chronic obstructive pulmonary disease|copd|"
              "inhalation|fibrosis|pneumonia|silicosis|asbestosis"),
}

CANDIDATES = '''
SELECT DISTINCT (CONCAT("AOP ", ?n) AS ?aopid) ?title WHERE {
  ?aop a aopo:AdverseOutcomePathway ; dc:title ?title ;
       aopo:has_adverse_outcome ?ao ; aopo:has_key_event ?ke .
  ?ao dc:title ?aotitle .
  BIND(CONCAT(LCASE(?title), " ", LCASE(?aotitle)) AS ?s)
  BIND(STRAFTER(STR(?aop), "aop/") AS ?n)
  FILTER regex(?s, "%s", "i") }
'''


def sparql(query: str) -> list[dict]:
    body = urllib.parse.urlencode({"query": query, "format": "json"}).encode()
    req = urllib.request.Request(
        m1.ENDPOINT, data=body,
        headers={"Accept": "application/sparql-results+json"})
    with urllib.request.urlopen(req, timeout=300) as fh:
        return json.load(fh)["results"]["bindings"]


def her_kes(organ: str) -> set[str]:
    path = NETWORKS / f"{organ}_nodes.tsv"
    rows = csv.DictReader(path.read_text(encoding="utf-8-sig").splitlines(), delimiter="\t")
    return {r["id"].strip() for r in rows}


def write_list(organ: str, aops: list[tuple[str, str]], provenance: str) -> None:
    path = NETWORKS / f"sidra_aops_{organ}.tsv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["aopid", "title"])
        for aopid, title in aops:
            w.writerow([aopid, title])
    print(f"  -> {path.relative_to(m1.ROOT)}  ({len(aops)} AOPs)  [{provenance}]")


def main() -> None:
    m1.verify_pin()

    membership: dict[str, set[str]] = {}
    titles: dict[str, str] = {}
    for b in sparql('SELECT DISTINCT (CONCAT("AOP ", ?n) AS ?aopid) ?keid ?title WHERE { '
                    '?aop a aopo:AdverseOutcomePathway ; dc:title ?title ; '
                    'aopo:has_key_event ?ke . ?ke rdfs:label ?keid . '
                    'BIND(STRAFTER(STR(?aop), "aop/") AS ?n) }'):
        aopid = b["aopid"]["value"]
        membership.setdefault(aopid, set()).add(b["keid"]["value"])
        titles[aopid] = b["title"]["value"]

    m1.banner("liver — from her own kes.tsv")
    rows = csv.DictReader((NETWORKS / "kes.tsv").read_text(encoding="utf-8-sig").splitlines(),
                          delimiter="\t")
    liver = sorted({r["aopNumber"].strip() for r in rows},
                   key=lambda a: int(a.split()[1]))
    write_list("liver", [(a, titles.get(a, "")) for a in liver], "kes.tsv")

    for organ, kw in KEYWORDS.items():
        m1.banner(f"{organ} — recovered from the keyword search")
        cand = sorted({b["aopid"]["value"] for b in sparql(CANDIDATES % kw)},
                      key=lambda a: int(a.split()[1]))
        hers = her_kes(organ)
        keep = [a for a in cand if membership.get(a, set()) <= hers]
        union: set[str] = set().union(*(membership[a] for a in keep)) if keep else set()

        print(f"  {len(cand)} candidates -> {len(keep)} kept, {len(cand) - len(keep)} excluded")
        if union != hers:
            raise SystemExit(
                f"ABORT: recovered {organ} AOPs give {len(union)} KEs but her node table has "
                f"{len(hers)}. Missing {sorted(hers - union)}, extra {sorted(union - hers)}. "
                "The AOP list cannot be recovered unambiguously — ask Sidra for it.")
        print(f"  check passed: KE union == her node table exactly ({len(union)} KEs)")
        write_list(organ, [(a, titles.get(a, "")) for a in keep],
                   "recovered, verified against her node table")


if __name__ == "__main__":
    main()
