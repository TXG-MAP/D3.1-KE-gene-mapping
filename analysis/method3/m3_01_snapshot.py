#!/usr/bin/env python3
"""Stage 1 — snapshot everything the run will be adjudicated against.

Three sources, all captured to disk so the rest of the pipeline is offline and reproducible:

  1. The pinned AOP-Wiki release  -> KE titles, biological levels, full description text.
  2. The molAOP Builder's corpora -> its KE snapshot and its ~803-pathway WikiPathways corpus.
  3. The builder's approved mappings -> what already exists, so nothing is proposed twice.

The builder has no release pin of its own; this
stage closes it by recording the pathway count, the KE-snapshot size and the GMT export's
`x-export-revision` header in the manifest. Quote those three whenever a Method 3 number is
quoted, the way `pav:createdOn` is quoted for Method 1.

Read-only: one SPARQL SELECT and four GETs.
"""
import csv
import hashlib
import io

import pandas as pd

import m3lib as L  # must precede m1lib: importing it puts analysis/method1 on sys.path
import m1lib  # noqa: E402


def fetch_ke_descriptions() -> pd.DataFrame:
    """KE title, biolevel and description text from the pinned release."""
    m1lib.verify_pin()
    df = m1lib.run_query(m1lib.read_query("method3a_ke_descriptions.rq"))
    df = df.rename(columns={"ketitle": "title"})
    df["description"] = df["description"].fillna("")
    df["biolevel"] = df["biolevel"].fillna("")
    # One row per KE. A KE with two biolevels would duplicate; assert it does not happen.
    dupes = df["keid"].duplicated().sum()
    if dupes:
        raise SystemExit(f"{dupes} duplicated KE rows from the endpoint — investigate before "
                         "using this as the biolevel source, it decides the +1.0 bonus")
    return df


def main():
    L.banner("stage 1: snapshot the release, the builder corpora and the live mappings")

    # --- 1. the pinned release ---------------------------------------------------------
    print("\n  AOP-Wiki RDF (pinned):")
    kes = fetch_ke_descriptions()
    described = (kes["description"].str.len() > 0).sum()
    print(f"    {len(kes)} KEs, {described} with a description, "
          f"{(kes['biolevel'] != '').sum()} with a biological level")
    L.save_table(kes, L.CORPUS_DIR / "ke_descriptions.csv", sep=",")

    # --- 2. the builder's corpora ------------------------------------------------------
    print(f"\n  molAOP Builder ({L.BASE}):")
    ke_options = L.get("/get_ke_options", pace=2.0)
    L.dump_json(ke_options, L.CORPUS_DIR / "ke_options.json")
    print(f"    /get_ke_options       {len(ke_options)} KEs")

    pathways = L.get("/get_pathway_options", pace=2.0)
    L.dump_json(pathways, L.CORPUS_DIR / "pathway_options.json")
    # The [10, 500] gene-set-size filter did not go away when #238 was fixed — it MOVED. It
    # used to be baked into the corpus, so an excluded pathway was invisible everywhere; it
    # now travels as a per-pathway `inSuggestionCorpus` flag applied only in the ranker. So a
    # pathway outside the bounds is fully selectable by hand and will never be suggested. That
    # distinction matters for the report: WP699 (Aflatoxin B1 metabolism, 7 genes) is no longer
    # "blocked", it is "invisible to the ranker", which is a different recommendation.
    rankable = sum(1 for p in pathways if p.get("inSuggestionCorpus", True))
    print(f"    /get_pathway_options  {len(pathways)} selectable, {rankable} rankable "
          f"({len(pathways) - rankable} outside [10, 500] genes: selectable, never suggested)")

    # --- 3. what is already approved ---------------------------------------------------
    mappings, page = [], 1
    while True:
        payload = L.get("/api/v1/mappings", {"page": page, "per_page": 200}, pace=1.0)
        data = payload.get("data") or []
        mappings.extend(data)
        if len(data) < 200:
            break
        page += 1
    L.dump_json({"data": mappings}, L.CORPUS_DIR / "live_mappings.json")
    print(f"    /api/v1/mappings      {len(mappings)} approved KE-WP pairs")

    body, headers = L.get("/exports/gmt/ke-wp", raw=True, pace=1.0)
    (L.CORPUS_DIR / "ke-wp.gmt").write_bytes(body)
    export_revision = headers.get("x-export-revision", "")
    gmt_lines = [line for line in body.decode().splitlines()
                 if line.strip() and not line.startswith("#")]
    gmt_sets = len(gmt_lines)

    # The file's own sha256 is NOT a usable pin: the export stamps a fresh `# generated:` line
    # and does not order genes within a set stably, so two exports of identical content hash
    # differently. Verified 2026-08-07 — re-exporting gave the same 246 sets, the same titles
    # and the same gene membership on every one, a different file sha256, and the SAME
    # x-export-revision. So the header is the content identifier, and this order-independent
    # digest is the local check that agrees with it. Keep both: if they ever disagree, the
    # header is being computed over something other than the payload.
    gmt_content_sha = hashlib.sha256("\n".join(sorted(
        "\t".join([f.split("\t")[0], f.split("\t")[1]] + sorted(f.split("\t")[2:]))
        for f in gmt_lines)).encode()).hexdigest()
    print(f"    /exports/gmt/ke-wp    {gmt_sets} sets, revision {export_revision}")

    # The GMT and the API must agree; if they do not, one of them is being read wrongly and
    # every downstream coverage number inherits the error.
    if gmt_sets != len(mappings):
        print(f"    ! GMT sets ({gmt_sets}) != API mappings ({len(mappings)}) — "
              "check the join before quoting coverage")

    # --- reconcile the builder's KE snapshot against the pin ---------------------------
    builder_kes = {k.get("KElabel") for k in ke_options}
    release_kes = set(kes["keid"])
    missing = sorted(release_kes - builder_kes, key=lambda x: int(x.split()[1]))
    print(f"\n  KE snapshot drift: the builder serves {len(builder_kes)} of the release's "
          f"{len(release_kes)} KEs; {len(missing)} are absent (upstream issue #239)")

    with open(L.CORPUS_DIR / "kes_absent_from_builder.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["keid", "title"])
        titles = dict(zip(kes["keid"], kes["title"]))
        for k in missing:
            w.writerow([k, titles.get(k, "")])
    print(f"  -> {(L.CORPUS_DIR / 'kes_absent_from_builder.tsv').relative_to(L.ROOT)}")

    L.write_manifest("01_snapshot", {
        "release_kes": len(release_kes),
        "release_kes_described": int(described),
        "builder_ke_snapshot": len(builder_kes),
        "builder_kes_absent_vs_release": len(missing),
        "builder_pathways_selectable": len(pathways),
        "builder_pathways_rankable": rankable,
        "builder_pathways_selectable_not_rankable": len(pathways) - rankable,
        "approved_wp_mappings": len(mappings),
        "gmt_sets": gmt_sets,
        "gmt_export_revision": export_revision,
        "gmt_sha256": L.sha256(L.CORPUS_DIR / "ke-wp.gmt"),
        "gmt_content_sha256": gmt_content_sha,
        "gmt_sha256_note": "gmt_sha256 is the raw file and is NOT stable across exports of "
                           "identical content (fresh timestamp, unstable gene order). Pin on "
                           "gmt_export_revision, or on gmt_content_sha256 which is computed "
                           "over sorted sets and sorted genes.",
        "note": "The builder has no upstream release pin. These five counts plus the export "
                "revision ARE the Method 3 corpus pin — quote them with any Method 3 number.",
    })


if __name__ == "__main__":
    main()
