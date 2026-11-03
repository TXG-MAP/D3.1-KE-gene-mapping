#!/usr/bin/env python3
"""Stage 12 — Method 3 coverage of the four frozen organ networks, on the closed queue.

This is the stage that turns "223 proposals lodged" into a reportable coverage number. It could
not run until the curators had worked the queue: a lodged proposal is not a mapping, and the
brain/kidney/lung figures land *below* the local run's selections wherever one was rejected.

**The queue closed on 2026-08-13** (Marvin): every one of the 223 proposals lodged on 2026-08-10
has been adjudicated by a curator. What is absent from the approved set is therefore a
rejection, not a backlog — which is what makes the accept rate below a real number rather than
a snapshot of how far Sidra had got.

Scope here is all **four** organ networks, not the three of stages 0-11. Liver was curated in
July 2026 (D3.2) and its coverage is re-derived here on the frozen frame, replacing the figure
computed against an earlier liver AOP selection — both frames hold exactly 196 KEs but
they are not the same 196, which is precisely how the mismatch went unnoticed (claim C09).

Two sources, deliberately both:

  * ``/api/v1/mappings`` — structured, paged, and the authority for which pairs are approved.
  * ``/exports/gmt/ke-wp`` — the same content plus gene symbols, and the artefact a reader can
    download. Parsed only to cross-check the pair count and to carry the gene sets to stage 13.

Coverage counts **every** approved pair on a network's KEs, not only the ones this campaign
lodged: 94 of liver's 246 predated the July campaign, and a reader asking "how much of the
liver network does Method 3 cover" wants the whole approved set.

Network: 2 GETs plus pagination. Writes nothing to the builder.
"""
from __future__ import annotations

import csv
import hashlib
from datetime import date

import pandas as pd

import m3lib as L

# Stage 12 covers all four networks; m3lib.ORGANS is the three this run curated.
FOUR = ["brain", "kidney", "liver", "lung"]
NODE_FILE = {"brain": "brain_nodes.tsv", "kidney": "kidney_nodes.tsv",
             "liver": "liver_nodes.tsv", "lung": "lungs_nodes.tsv"}
EXPECTED_KES = {"brain": 203, "kidney": 68, "liver": 196, "lung": 134}

#: The union of the four frozen networks. 538 distinct KEs, 50 of which sit in more than one
#: network — so per-organ rows must never be summed.
EXPECTED_UNION = 538

LEDGER = L.M3 / "submission" / "ledger.tsv"
STAMP = date.today().isoformat()


def load_frames() -> dict[str, pd.DataFrame]:
    """The four frozen organ networks, as node tables keyed by KE id."""
    frames = {}
    for organ in FOUR:
        path = L.ROOT / "data" / "networks" / NODE_FILE[organ]
        df = pd.read_csv(path, sep="\t")
        df["id"] = df["id"].str.strip()
        if len(df) != EXPECTED_KES[organ]:
            raise SystemExit(
                f"{organ}: {len(df)} KEs, expected {EXPECTED_KES[organ]}. The frozen networks "
                f"have moved; Q9 was settled on these files and every method reports over them."
            )
        frames[organ] = df
    union = set().union(*(set(df["id"]) for df in frames.values()))
    if len(union) != EXPECTED_UNION:
        raise SystemExit(f"union is {len(union)} KEs, expected {EXPECTED_UNION}")
    # 50 KEs sit in more than one network but they account for 63 excess rows, because 10 of
    # them sit in three or more. Quote the 50; the 63 is not a count of anything reportable.
    counts = pd.Series([i for df in frames.values() for i in df["id"]]).value_counts()
    shared = int((counts > 1).sum())
    if shared != 50:
        raise SystemExit(f"{shared} KEs in more than one network, expected 50 (Q9)")
    print(f"  four frozen networks: {sum(len(d) for d in frames.values())} rows, "
          f"{len(union)} distinct KEs, {shared} of them in more than one network")
    return frames


def fetch_mappings() -> list[dict]:
    """Every approved KE-WP pair the service holds, paged."""
    mappings, page = [], 1
    while True:
        payload = L.get("/api/v1/mappings", {"page": page, "per_page": 200}, pace=1.0)
        data = payload.get("data") or []
        mappings.extend(data)
        if len(data) < 200:
            break
        page += 1
    return mappings


def fetch_gmt() -> tuple[list[str], str, str]:
    """The GMT export. Returns (data lines, export revision, order-independent content sha).

    The file's own sha256 is not a pin — the export stamps a fresh `# generated:` line and does
    not order genes stably within a set. The header revision is the content identifier; this
    digest is the local check that should agree with it.
    """
    body, headers = L.get("/exports/gmt/ke-wp", raw=True, pace=1.0)
    (L.CORPUS_DIR / f"ke-wp_final_{STAMP}.gmt").write_bytes(body)
    lines = [ln for ln in body.decode().splitlines() if ln.strip() and not ln.startswith("#")]
    content_sha = hashlib.sha256("\n".join(sorted(
        "\t".join([f.split("\t")[0], f.split("\t")[1]] + sorted(f.split("\t")[2:]))
        for f in lines)).encode()).hexdigest()
    return lines, headers.get("x-export-revision", ""), content_sha


def ledger_outcomes(approved_pairs: set[tuple[str, str]]) -> pd.DataFrame:
    """What became of the 223 proposals lodged on 2026-08-10.

    With the queue closed, absence from the approved set is a rejection. While it was open this
    same computation could only say "not yet approved", which is why nothing quoted it.
    """
    rows = []
    with open(LEDGER, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ke, wp = r["ke_id"].strip(), r["wp_id"].strip()
            rows.append({
                "proposal_id": int(r["proposal_id"]),
                "ke_id": ke,
                "wp_id": wp,
                "predicted_confidence": r["confidence"],
                "step1": r["step1"], "step2": r["step2"],
                "step3": r["step3"], "step4": r["step4"],
                "outcome": "approved" if (ke, wp) in approved_pairs else "rejected",
            })
    return pd.DataFrame(rows).sort_values("proposal_id")


def coverage(frames: dict[str, pd.DataFrame],
             pairs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-organ coverage, plus the per-KE long table it is computed from."""
    per_ke, rows = [], []
    for organ in FOUR:
        frame = frames[organ]
        sub = pairs[pairs["ke_id"].isin(set(frame["id"]))]
        by_ke = sub.groupby("ke_id")["wp_id"].nunique()
        mapped = len(by_ke)
        rows.append({
            "organ": organ,
            "kes_in_frame": len(frame),
            "kes_with_approved_mapping": mapped,
            "pct_kes_covered": round(100 * mapped / len(frame), 1),
            "kes_without": len(frame) - mapped,
            "approved_pairs": len(sub),
            "distinct_pathways": sub["wp_id"].nunique(),
            "median_pathways_per_mapped_ke": float(by_ke.median()) if mapped else 0.0,
        })
        for _, ke in frame.iterrows():
            hits = sub[sub["ke_id"] == ke["id"]]
            per_ke.append({
                "organ": organ,
                "ke_id": ke["id"],
                "ke_name": ke["name"],
                "ke_type": ke["type"],
                "n_approved_pathways": len(hits),
                "pathways": "|".join(sorted(hits["wp_id"])),
            })

    # The union row is computed over distinct KEs, never by summing the four above: 50 KEs sit
    # in more than one network and would otherwise be counted twice.
    union_kes = set().union(*(set(df["id"]) for df in frames.values()))
    usub = pairs[pairs["ke_id"].isin(union_kes)]
    umapped = usub["ke_id"].nunique()
    rows.append({
        "organ": "four-organ union",
        "kes_in_frame": len(union_kes),
        "kes_with_approved_mapping": umapped,
        "pct_kes_covered": round(100 * umapped / len(union_kes), 1),
        "kes_without": len(union_kes) - umapped,
        "approved_pairs": len(usub),
        "distinct_pathways": usub["wp_id"].nunique(),
        "median_pathways_per_mapped_ke": float(
            usub.groupby("ke_id")["wp_id"].nunique().median()) if umapped else 0.0,
    })
    return pd.DataFrame(rows), pd.DataFrame(per_ke)


def markdown(cov: pd.DataFrame, outcomes: pd.DataFrame, meta: dict) -> str:
    out = [
        "# Method 3 — coverage of the four frozen organ networks",
        "",
        f"Generated {meta['generated']} from the molAOP Builder's approved mappings "
        f"(export revision `{meta['revision']}`, {meta['pairs']} approved KE-WP pairs over "
        f"{meta['kes']} Key Events service-wide).",
        "",
        "Frames are the frozen networks in `data/networks/*_nodes.tsv` (Q9). **The four organ "
        "rows must not be summed** — 50 KEs sit in more than one network; the union row already "
        "accounts for them.",
        "",
        "| Organ | KEs | with a mapping | coverage | without | approved pairs | distinct WPs |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, r in cov.iterrows():
        bold = "**" if r["organ"] == "four-organ union" else ""
        out.append(
            f"| {bold}{r['organ'].capitalize() if not bold else r['organ']}{bold} "
            f"| {r['kes_in_frame']} | {bold}{r['kes_with_approved_mapping']}{bold} "
            f"| {bold}{r['pct_kes_covered']} %{bold} | {r['kes_without']} "
            f"| {r['approved_pairs']} | {r['distinct_pathways']} |"
        )

    n = len(outcomes)
    app = int((outcomes["outcome"] == "approved").sum())
    rej = n - app
    out += [
        "",
        "## Outcome of the 223 proposals lodged on 2026-08-10",
        "",
        f"The queue closed on {meta['queue_closed']}: every proposal has been adjudicated, so "
        f"absence from the approved set is a rejection rather than a backlog.",
        "",
        f"| | |", "| --- | ---: |",
        f"| Lodged (ids 333-555) | {n} |",
        f"| Approved | **{app}** |",
        f"| Rejected | **{rej}** |",
        f"| Accept rate | **{100 * app / n:.1f} %** |",
        "",
        "By the confidence tier predicted offline before submission:",
        "",
        "| Predicted tier | lodged | approved | rejected | accept rate |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for tier in ("high", "medium", "low"):
        sub = outcomes[outcomes["predicted_confidence"] == tier]
        if not len(sub):
            continue
        a = int((sub["outcome"] == "approved").sum())
        out.append(f"| {tier} | {len(sub)} | {a} | {len(sub) - a} "
                   f"| {100 * a / len(sub):.1f} % |")

    rejected = outcomes[outcomes["outcome"] == "rejected"]
    if len(rejected):
        out += ["", "Rejected proposals:", "",
                "| id | KE | WP | predicted tier |", "| ---: | --- | --- | --- |"]
        for _, r in rejected.iterrows():
            out.append(f"| {r['proposal_id']} | {r['ke_id']} | {r['wp_id']} "
                       f"| {r['predicted_confidence']} |")
    return "\n".join(out) + "\n"


def main():
    L.banner("stage 12: Method 3 coverage on the four frozen networks")

    frames = load_frames()

    print("\n  fetching the approved set")
    mappings = fetch_mappings()
    L.dump_json({"data": mappings}, L.CORPUS_DIR / f"live_mappings_final_{STAMP}.json")
    pairs = pd.DataFrame([{
        "ke_id": m["ke_id"].strip(),
        "ke_name": m.get("ke_name", ""),
        "ke_bio_level": m.get("ke_bio_level", ""),
        "wp_id": m["pathway_id"].strip(),
        "wp_title": m.get("pathway_title", ""),
        "confidence": m.get("confidence_level", ""),
        "approved_at": (m.get("provenance") or {}).get("approved_at", ""),
    } for m in mappings])
    print(f"    /api/v1/mappings      {len(pairs)} approved pairs over "
          f"{pairs['ke_id'].nunique()} KEs")

    gmt_lines, revision, content_sha = fetch_gmt()
    print(f"    /exports/gmt/ke-wp    {len(gmt_lines)} sets, revision {revision}")
    if len(gmt_lines) != len(pairs):
        raise SystemExit(
            f"GMT holds {len(gmt_lines)} sets against {len(pairs)} API mappings. One of the two "
            f"is being read wrongly and every coverage number below would inherit the error."
        )

    dup = pairs.duplicated(subset=["ke_id", "wp_id"]).sum()
    if dup:
        raise SystemExit(f"{dup} duplicate KE-WP pairs in the approved set")

    print("\n  outcome of the lodged proposals")
    if LEDGER.exists():
        outcomes = ledger_outcomes(set(zip(pairs["ke_id"], pairs["wp_id"])))
    else:
        # The submission ledger is not part of the public repository; the outcome table it
        # produced is, and is reused as-is.
        outcomes = pd.read_csv(L.DERIVED_DIR / "method3_proposal_outcomes.tsv", sep="\t")
    app = int((outcomes["outcome"] == "approved").sum())
    print(f"    {len(outcomes)} lodged -> {app} approved, {len(outcomes) - app} rejected "
          f"({100 * app / len(outcomes):.1f} % accepted)")

    print("\n  coverage")
    cov, per_ke = coverage(frames, pairs)
    for _, r in cov.iterrows():
        print(f"    {r['organ']:<18} {r['kes_with_approved_mapping']:>4}/{r['kes_in_frame']:<4} "
              f"({r['pct_kes_covered']:>5.1f} %)  {r['approved_pairs']:>4} pairs")

    L.save_table(cov, L.DERIVED_DIR / "method3_coverage_by_organ.tsv")
    L.save_table(per_ke, L.DERIVED_DIR / "method3_coverage_by_ke.tsv")
    L.save_table(outcomes, L.DERIVED_DIR / "method3_proposal_outcomes.tsv")
    L.save_table(pairs, L.DERIVED_DIR / "method3_approved_pairs.tsv")

    meta = {
        "generated": STAMP,
        "revision": revision,
        "pairs": len(pairs),
        "kes": int(pairs["ke_id"].nunique()),
        "queue_closed": "2026-08-13",
    }
    md = L.ROOT / "tables" / "method3_coverage.md"
    md.parent.mkdir(parents=True, exist_ok=True)
    md.write_text(markdown(cov, outcomes, meta))
    print(f"  -> {md.relative_to(L.ROOT)}")

    L.write_manifest("12_coverage", {
        "purpose": "Method 3 coverage of the four frozen organ networks, on the closed queue",
        "queue_closed_on": "2026-08-13",
        "queue_closed_by": "curator adjudication complete (Marvin, 2026-08-13)",
        "gmt_export_revision": revision,
        "gmt_content_sha256": content_sha,
        "approved_pairs": len(pairs),
        "approved_kes": int(pairs["ke_id"].nunique()),
        "frames": {o: int(len(frames[o])) for o in FOUR},
        "frame_union_kes": EXPECTED_UNION,
        "proposals_lodged": len(outcomes),
        "proposals_approved": app,
        "proposals_rejected": len(outcomes) - app,
        "coverage": cov.to_dict("records"),
        "note": "Coverage counts every approved pair on a network's KEs, including mappings "
                "that predate the D3.1/D3.2 campaigns. Organ rows must not be summed: 50 KEs "
                "sit in more than one frozen network.",
        "outputs": [
            "data/method3/derived/method3_coverage_by_organ.tsv",
            "data/method3/derived/method3_coverage_by_ke.tsv",
            "data/method3/derived/method3_proposal_outcomes.tsv",
            "data/method3/derived/method3_approved_pairs.tsv",
            f"data/method3/corpus/ke-wp_final_{STAMP}.gmt",
            f"data/method3/corpus/live_mappings_final_{STAMP}.json",
            "tables/method3_coverage.md",
        ],
    })


if __name__ == "__main__":
    main()
