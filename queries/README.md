# SPARQL queries

Every query used in the report lives here as a runnable `.rq` file and is reproduced in
Appendix A of the report. The files and the appendix are kept in sync.

## Conventions

- One query per file, named `<method>_<what-it-returns>.rq`.
- A header comment block gives: endpoint, RDF release the results in the report came from,
  which script or notebook calls it, and any non-obvious modelling choice (e.g. why a block is
  `OPTIONAL`).
- Placeholders that the calling code substitutes are written as a quoted literal with a
  trailing comment, so the file stays syntactically valid and can be pasted into SNORQL as-is.

## Endpoints

Two endpoints, split by purpose:

| Endpoint | Use for | Pinned to |
| --- | --- | --- |
| <https://aopwiki.rdf.bigcat-bioinformatics.org/sparql/> | **every content number** — KE inventory, all gene mappings (Method 1), KE-component counts | `97191db`, generated 2026-08-01 |
| <https://aopwiki-multirdf.vhp4safety.nl/sparql> | **temporal trends only** — Method 2 Figure 2 C/D/E | 34 quarterly graphs, `2018-04-01` … `2026-07-01` |

**Method 1 must never source from the multi-version endpoint.** It is built by a different
converter that skips the gene-mapping step — no BERN2 enrichment, 382 KEs with genes against
488. Both pins are machine-readable in
`pinned-release.yaml`.

The endpoint regenerates weekly and will drift off the pinned commit. Check with:

```sparql
PREFIX pav: <http://purl.org/pav/>
SELECT ?s ?created WHERE { ?s pav:createdOn ?created }   # expect 2026-08-01
```

To distinguish two runs from the same date, match the `void:triples` counts —
138968 / 8360 / 214561.

## Inventory

| File | Purpose | Status |
| --- | --- | --- |
| `method1_ke_gene_annotations.rq` | KE → NCBI Gene per AOP (union of regex + NER) | superseded by `method1b`; kept because it is the query reproduced in the draft |
| `method1b_regex_vs_ner.rq` | KE → gene with per-gene detector provenance | **in use** — `m1_01_fetch.py` |
| `method1c_ke_inventory.rq` | every KE with title and biological level (the denominator) | **in use** — `m1_01_fetch.py` |
| `method1d_ker_gene_annotations.rq` | KER → genes, the layer the report does not count | **in use** — `m1_01_fetch.py` |
| `method1e_gene_xrefs.rq` | HGNC → Ensembl / NCBI Gene via `owl:sameAs` | **in use** — `m1_03_harmonise.py` |
| `method2_kec_ontology_sources.rq` | Ontology source distribution for BP/BO components | to write |
| `method2_ke_biological_process.rq` | KE → GO BP term (Process slot) | to write |
| `network_organ_aops.rq` | AOP → KE membership with role (MIE/KE/AO); the organ frame is applied downstream | **in use** — `m1_01_fetch.py` |
