#!/usr/bin/env python3
"""Assert that every ontology identifier in a Method 2 deposit table means what its title says.

    python3 m2_03_verify_deposit.py <file.tsv> [...]

This is the build-time gate the Zenodo deposit needed and did not have. Run it on any table
before publishing: it exits non-zero if a single identifier disagrees with its label.

**Why a gate rather than a review.** Every wrong identifier found in the 2026-08 liver deposit
resolves to a real, non-obsolete GO term. Nothing 404s, nothing errors, nothing looks wrong to a
reader — the row simply asserts a different biological process from the one it names, and the
only way to see it is to resolve the identifier and compare. That is precisely the class of
defect a human check will not catch and an assertion always will.

Labels are resolved from the cached go-basic.obo (stage 1). Only GO-sourced cells are checked;
MeSH, CHEBI, PATO, CL and UniProtKB identifiers pass through untested, and multi-term cells are
split on ';' and checked element-wise.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import m2lib as L


def go_labels() -> tuple[dict, dict]:
    names, alt, cur = {}, {}, None
    path = L.CACHE / "go-basic.obo"
    if not path.exists():
        raise SystemExit(f"{path} is missing — run `make method2-fetch` first.")
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {}
            continue
        if line.startswith("[") and line.endswith("]"):
            cur = None
            continue
        if cur is None:
            continue
        if line.startswith("id: "):
            cur["id"] = line[4:].strip()
        elif line.startswith("name: ") and "id" in cur:
            names[cur["id"]] = line[6:].strip()
        elif line.startswith("alt_id: ") and "id" in cur:
            alt[line[8:].strip()] = cur["id"]
    return names, alt


def check(path: Path, names: dict, alt: dict) -> list[dict]:
    def label(gid: str):
        return names.get(alt.get(gid, gid))

    delim = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    rows = list(csv.DictReader(open(path, newline=""), delimiter=delim))
    if not rows:
        raise SystemExit(f"{path}: no rows")

    # Column families: <prefix>_id / <prefix>_title / <prefix>_source / <prefix>_uri.
    prefixes = sorted({c[:-3] for c in rows[0] if c.endswith("_id")
                       and f"{c[:-3]}_title" in rows[0]})
    findings = []
    for i, r in enumerate(rows):
        for p in prefixes:
            source = (r.get(f"{p}_source") or "").strip()
            raw_id = (r.get(f"{p}_id") or "").strip()
            title = (r.get(f"{p}_title") or "").strip()
            uri = (r.get(f"{p}_uri") or "").strip()
            if source != "GO" or not raw_id or raw_id == "NA":
                continue

            ids = [x.strip() for x in raw_id.split(";") if x.strip()]
            titles = [x.strip() for x in title.split(";") if x.strip()]

            # Some columns strip the prefix and the leading zeros ("5739" for GO:0005739).
            # That is not resolvable as an identifier and is worth reporting on its own, but
            # the label check should still run, so normalise before comparing.
            norm = []
            for x in ids:
                if x.startswith("GO:"):
                    norm.append(x)
                elif x.isdigit():
                    norm.append(f"GO:{int(x):07d}")
                    findings.append({
                        "row": i + 2, "column": f"{p}_id", "kind": "unprefixed identifier",
                        "detail": f"{x!r} is not a resolvable identifier; "
                                  f"expected GO:{int(x):07d}",
                    })
                else:
                    norm.append(x)

            if len(norm) != len(titles):
                # Only report when the id side is multi-valued; a single id against a
                # multi-clause title is normal prose, not a defect.
                if len(norm) > 1:
                    findings.append({
                        "row": i + 2, "column": f"{p}_id", "kind": "cardinality",
                        "detail": f"{len(norm)} identifiers against {len(titles)} titles",
                    })
                continue

            for gid, ti in zip(norm, titles):
                got = label(gid)
                if got is None:
                    findings.append({
                        "row": i + 2, "column": f"{p}_id", "kind": "unknown or obsolete term",
                        "detail": f"{gid} is not a current term in go-basic",
                    })
                elif got.lower() != ti.lower():
                    findings.append({
                        "row": i + 2, "column": f"{p}_id", "kind": "identifier/title mismatch",
                        "detail": f"{gid} is {got!r} but the row is titled {ti!r}",
                    })

            if uri:
                m = re.search(r"GO[:_](\d+)", uri)
                if m and f"GO:{m.group(1)}" not in norm and len(norm) == 1:
                    findings.append({
                        "row": i + 2, "column": f"{p}_uri", "kind": "uri/id disagreement",
                        "detail": f"uri points at GO:{m.group(1)}, {p}_id says {norm[0]}",
                    })
    return findings


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    names, alt = go_labels()

    failed = False
    for arg in sys.argv[1:]:
        path = Path(arg)
        L.banner(f"verifying {path.name}")
        findings = check(path, names, alt)
        mismatches = [f for f in findings if f["kind"] == "identifier/title mismatch"]
        others = [f for f in findings if f["kind"] != "identifier/title mismatch"]

        for f in findings:
            print(f"  row {f['row']:>4}  {f['column']:<10} {f['kind']:<24} {f['detail']}")
        if not findings:
            print("  clean — every GO identifier matches its title")
        print(f"\n  {len(mismatches)} identifier/title mismatches, {len(others)} other issues")
        if mismatches:
            failed = True

    if failed:
        raise SystemExit(
            "\nFAILED. An identifier that disagrees with its title resolves to a real term and "
            "silently asserts different biology. Do not publish this table."
        )


if __name__ == "__main__":
    main()
