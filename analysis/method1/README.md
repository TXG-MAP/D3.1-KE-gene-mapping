# Method 1 — text-based KE→gene mapping, analysis pipeline

Five scripts, run in order. `make method1` from the repository root runs all of them;
`make method1-analyse` re-runs everything except the fetch, which is the one you want while
iterating on statistics or figures.

| Stage | Script | Talks to the network? | Writes |
| --- | --- | --- | --- |
| 1 | `m1_01_fetch.py` | **yes** — SPARQL | `data/method1/raw/*.csv` |
| 2 | `m1_02_frames.py` | no | `data/method1/frames/*.csv` |
| 3 | `m1_03_harmonise.py` | **yes** — BridgeDb (skippable) | `data/method1/derived/gene_id_harmonisation.csv` |
| 4 | `m1_04_stats.py` | no | `tables/method1_statistics.md`, `…/method1_statistics.json` |
| 5 | `m1_05_figures.py` | no | `figures/final/FIG07,08,09,12,13*` |

`m1lib.py` holds the shared configuration, the SPARQL client, and the plot style. Each stage
writes a `data/method1/manifest_<stage>.json` recording the release, the repository revision,
row counts and SHA-256 digests.

The mapping itself is **not** performed here. It happens upstream in the AOP-Wiki RDF
conversion pipeline (https://github.com/marvinm2/AOPWikiRDF); these scripts query
the result and analyse it.

## The three things that will bite

**The endpoint is rolling.** It regenerates weekly and will drift off the pinned commit
`97191db`. Stage 1 verifies `pav:createdOn` and the three `void:triples` counts before
fetching, and aborts on a mismatch. Overriding with `ALLOW_DRIFT=1` produces numbers that do
not correspond to the stated pin — re-pin `pinned-release.yaml` instead, or load the TTLs by
SHA from the `data_dumps` URLs.

**Booleans are not `true`.** Virtuoso serialises `xsd:boolean` as `1`/`0`. Parsing the
detector-provenance flags by matching the string `"true"` yields *every* flag `False`, which
reads as "BERN2 contributed nothing" rather than as an error. Stage 1 raises if no pair
carries a flag, because that pattern is a parse failure and never a data state.

**BridgeDb's batch endpoint ignores `dataSource`.** `/xrefs/` honours it and returns bare
identifiers; `/xrefsBatch/` ignores it and returns every system, prefixed (`En:ENSG…`,
`L:5468`, `Pd:`, `T:GO:…`). Filtering must be done client-side on the prefix, and `L:` has to
be matched exactly or it also swallows `Il:` (Illumina probes). Getting this wrong reports a
~100 % mapping rate for any target. Stage 3 also refuses to write a result below an 80 %
mapping rate, which is the guard the July 2026 incident argues for — a TSV→JSON change then
cut gene counts by 78.7 % on an HTTP 200.

## The frames, and what is blocked

An organ AOP network is a hand-made selection of AOPs, not something the RDF knows. Stage 2
materialises two:

- **`legacy2024`** — the February 2024 AO→organ classification, all four organs. Provisional;
  it is what the numbers currently in the draft were computed on, and it covers 360 of the
  release's 595 AOPs.
- **`sidra`** — the authoritative selection strategy (MS3.1).
  **Liver only.** Brain, kidney and lung have not been produced and this pipeline does not
  invent them; it records them as missing and every four-organ result is therefore on the
  provisional frame.

Running both is the point. It separates the **detector** effect (regex vs BERN2, one frame)
from the **frame** effect (two selections, one detector), which are confounded in any
comparison against the draft.

## Reproducing the draft's numbers

`m1_04_stats.py` block D compares against the 2025-02-19 figures carried in
the report's Method 1 section. Those are hard-coded in `DRAFT_NUMBERS` so the
comparison survives the section being rewritten — update them there if the draft changes.

## Regenerating a figure only

```bash
cd analysis/method1 && python3 m1_05_figures.py
```

Figures read `data/method1/derived/method1_statistics.json`, so restyling never needs the
endpoint. Colours live in `m1lib.py` and were checked for colour-vision-deficiency
separation; re-check them if they change.
