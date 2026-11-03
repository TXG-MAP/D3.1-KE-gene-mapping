# D3.1 — Evaluation and comparison of KE-gene mapping sets

Analysis code, queries, data extracts and figures behind deliverable D3.1 of the EFSA-funded
project TXG-MAP (GP/EFSA/ED/2022/01). The deliverable compares three ways of assigning gene sets
to the Key Events (KEs) of four organ-specific Adverse Outcome Pathway (AOP) networks (brain,
kidney, liver and lung):

| | Method 1 | Method 2 | Method 3 |
| --- | --- | --- | --- |
| Route | text mining of KE descriptions | Gene Ontology Biological Process annotation | KE → WikiPathways molecular pathways |
| Source | AOP-Wiki RDF (HGNC dictionary + BERN2) | curated KE–GO mappings, served by the molAOP Builder | approved KE–pathway mappings, served by the molAOP Builder |

All three are harmonised to Ensembl gene identifiers for the comparison.

## Pinned inputs

Every content number comes from one release of the AOP-Wiki RDF: commit `97191db` of
[marvinm2/AOPWikiRDF](https://github.com/marvinm2/AOPWikiRDF), data generated on 1 August 2026
(595 AOPs, 1 595 KEs, 2 359 KERs). The endpoint rolls weekly, so the pin and the immutable
Turtle URLs are recorded in [`pinned-release.yaml`](pinned-release.yaml), and the extract the
Method 1 numbers were computed from is kept in `data/method1/raw/`.

The molAOP Builder has no release pin of its own. Its KE snapshot, pathway corpus and approved
mapping exports, as read on 11–14 August 2026, are kept in `data/method3/corpus/` and
`data/method2/corpus/`.

## Layout

```
queries/        SPARQL queries, one per file (Appendix A of the report)
data/networks/  the four organ AOP networks (AOP lists, nodes, edges, KE inventory)
data/method1/   Method 1 extract, organ frames, derived tables and run manifests
data/method2/   Method 2 gene sets derived from the builder's KE–GO export
data/method3/   Method 3 builder snapshot, KE frame, coverage and gene-set tables
data/comparison/ cross-method tables
analysis/       pipeline scripts, by method; analysis/comparison/ for the cross-method stage
tables/         generated summary tables
figures/final/  figures with a provenance strip under the plot
figures/clean/  the same figures without it
```

## Running

Python 3.11 or later; `pip install -r requirements.txt`. Then, for example:

```bash
make method1-analyse      # Method 1 from the cached extract, no endpoint needed
make method3-coverage     # Method 3 coverage of the four networks
make builder-genesets     # Methods 2 and 3 gene sets from the builder exports
make comparison           # the cross-method comparison
```

`make help` lists every stage. Stages that call BridgeDb or a live endpoint accept
`OFFLINE=1` where an offline route exists. After any re-run against a live service, compare
the counts with the committed tables before trusting them.

`make method3-coverage` reads the molAOP Builder's live API, which keeps growing as curation
continues. The committed Method 3 tables record the state the report uses: the approved set
when the proposal queue closed on 13 August 2026 (432 pairs over 314 KEs of the four
networks). A re-run today will report more.

## Not in this repository

- **Method 2's curated KE–GO mapping tables** are published on Zenodo:
  [10.5281/zenodo.21264952](https://doi.org/10.5281/zenodo.21264952) (concept DOI, resolves to
  the latest version).
- **The Method 3 curation itself**, meaning the review of candidate pathways and their
  submission to the molAOP Builder. The analysis here starts from the builder's export of
  approved mappings (<https://molaop-builder.vhp4safety.nl/exports/gmt/ke-wp>), snapshotted in
  `data/method3/corpus/`.
- **The report text**, which is delivered to EFSA separately.

## Licence

Code is released under the MIT licence ([`LICENSE`](LICENSE)). Data tables and figures are
released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). AOP-Wiki content is
subject to the AOP-Wiki terms, and WikiPathways content to CC0.

## Citation

See [`CITATION.cff`](CITATION.cff).
