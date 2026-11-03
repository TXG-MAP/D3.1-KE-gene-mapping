# Method 1 — statistics

Generated from AOP-Wiki RDF commit `97191db`, data generated 2026-08-01, via `analysis/method1/m1_04_stats.py`.

Regenerate with `make method1`. Do not hand-edit: this file is overwritten.


## A. Repository-wide, all 1595 Key Events

Method 1 annotates at least one gene onto **488 of 1595 Key Events (30.6%, 95 % CI 28.4%–32.9%)**, using 1914 distinct genes over 4366 KE-gene pairs.


### A0. The denominator: which Key Events can this method reach at all?

Method 1 finds gene names by reading a Key Event's free-text description. A Key Event with no description cannot receive a gene annotation from it — that is not a miss, it is out of scope. So `30.6 % of Key Events` measures AOP-Wiki's description completeness at least as much as it measures the method, and the two should be separated.

| population | KEs | annotated | coverage % | CI low | CI high |
| --- | ---: | ---: | ---: | ---: | ---: |
| all Key Events in the release | 1595 | 488 | 30.6 | 28.4 | 32.9 |
| Key Events carrying a description | 635 | 488 | 76.9 | 73.4 | 80.0 |
| Key Events with no description | 960 | 0 | 0.0 |  |  |

**All 488 gene-annotated Key Events carry a description, and none of the 960 undescribed Key Events carries a gene.** The containment is exact, which confirms the description field is the sole route by which a KE acquires genes under this method and is asserted on every pipeline run.

So the method's own hit rate is **76.9% (488 of 635, 95 % CI 73.4%–80.0%)**, not 30.6 %. The two numbers answer different questions and the report needs both: 77% is how well the text mining works, 31% is how much of the AOP-Wiki it can currently serve. The gap between them is the 960 Key Events (60%) that have no description to mine — a gap that no improvement to the method can close, and a concrete curation recommendation for AOP-Wiki rather than a limitation of Method 1.


### A0b. How much text a Key Event has predicts what the method finds

| description length | KEs | annotated | median_length | median_genes | % annotated |
| --- | ---: | ---: | ---: | ---: | ---: |
| Q1 shortest | 160 | 66 | 404.5 | 0.0 | 41.2 |
| Q2 | 158 | 129 | 1067.5 | 3.0 | 81.6 |
| Q3 | 158 | 141 | 2078.0 | 6.0 | 89.2 |
| Q4 longest | 159 | 152 | 4616.0 | 13.0 | 95.6 |

Among described Key Events, description length correlates with the number of genes found at **Spearman rho = 0.67** (p = 6.3e-84). Annotated KEs have a median description of 1877 characters against 540 for unannotated ones (Mann-Whitney p = 9.7e-33). In the longest quartile the method annotates 96 % of Key Events; in the shortest, 41 %.

This is worth stating plainly because it constrains how the output should be read: the size of a Key Event's Method 1 gene set reflects how much has been written about it as much as how many genes are biologically involved. A KE with 13 genes and one with 3 are not necessarily different in biology — they may differ in curation effort. That makes raw gene-set size a poor basis for comparing Key Events, and it is a caution the comparison in Section 07 has to carry.


### A1. What BERN2 changed: the detector decomposition

The same comparison counted three ways. These are different questions and the report should say which one it is quoting: a gene can be found by both detectors on one KE and by only one on another, so the gene-level and pair-level splits do not agree and neither is wrong.

| counted over | regex (HGNC) | BERN2 (NER) | regex only | both | BERN2 only | union |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Key Events | 384 | 462 | 26 | 358 | 104 | 488 |
| unique genes | 901 | 1652 | 262 | 639 | 1013 | 1914 |
| KE-gene pairs | 1857 | 3670 | 696 | 1161 | 2509 | 4366 |


### A2. Coverage under each detector alone

What the deliverable would have reported had only one detector run. The regex column is the closest available reconstruction of the pre-2026-06 pipeline on the pinned release — it is *not* the number in the current draft, which was computed on an older AOP-Wiki release as well as an older detector.

| detector | KEs with >=1 gene | % of 1595 KEs | unique genes | KE-gene pairs |
| --- | ---: | ---: | ---: | ---: |
| regex / HGNC dictionary only | 384 | 24.08 | 901 | 1857 |
| BERN2 NER only | 462 | 28.97 | 1652 | 3670 |
| union (as released) | 488 | 30.60 | 1914 | 4366 |

Adding BERN2 to the dictionary matcher increased annotated Key Events by **+27.1%** and unique genes by **+112.4%**. BERN2 alone recovers 70.9% of what the dictionary found, so the two are complementary rather than one superseding the other: 262 genes are found *only* by the dictionary.

**McNemar exact test** on KE-level coverage (paired over all 1595 KEs): 104 KEs annotated only by BERN2 against 26 only by the dictionary, p = 3.03e-12. The detectors are not interchangeable.

**Wilcoxon signed-rank** on genes per KE, BERN2 against dictionary, over the 488 annotated KEs: median 5 vs 2 genes, W = 7044, p = 1.34e-53.


### A3. Genes per Key Event

| detector | n | total | mean | sd | median | q1 | q3 | min | max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| union (as released) | 488 | 4366.0 | 8.9 | 9.9 | 6.0 | 2.0 | 11.2 | 1.0 | 60.0 |
| BERN2 | 488 | 3670.0 | 7.5 | 8.5 | 5.0 | 2.0 | 10.0 | 0.0 | 56.0 |
| regex | 488 | 1857.0 | 3.8 | 5.0 | 2.0 | 1.0 | 5.0 | 0.0 | 29.0 |

Denominator is the 488 KEs carrying at least one gene, so the regex row counts zero for KEs only BERN2 found. The distribution is heavily right-skewed — the mean sits well above the median and the largest single KE carries 60 genes — which is why medians are reported throughout and why a mean gene count per organ would be misleading.


### A4. Coverage by biological level of organisation

Reported against both denominators, because the obvious objection to the raw gradient is that molecular Key Events might simply be better described. They are not: the gradient is a property of the method, and conditioning on having a description makes it stronger rather than weaker.

| biological level | KEs | described | annotated | % described | % annotated of all | % annotated of described | median genes / annotated KE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Molecular | 448 | 191 | 171 | 42.6 | 38.2 | 89.5 | 6.0 |
| Cellular | 539 | 213 | 175 | 39.5 | 32.5 | 82.2 | 7.0 |
| Tissue | 221 | 85 | 63 | 38.5 | 28.5 | 74.1 | 6.0 |
| Organ | 195 | 88 | 49 | 45.1 | 25.1 | 55.7 | 2.0 |
| Individual | 156 | 48 | 25 | 30.8 | 16.0 | 52.1 | 3.0 |
| Population | 36 | 10 | 5 | 27.8 | 13.9 | 50.0 | 1.0 |

Three chi-square tests of independence across the six levels:

- **coverage over all KEs**: chi2(5) = 36.5, p = 7.48e-07
- **whether a KE is described at all**: chi2(5) = 11.5, p = 0.0427
- **coverage among described KEs**: chi2(5) = 63.8, p = 2.03e-12

Description availability does vary by level, but weakly and **non-monotonically** — Organ Key Events are the most often described (45 %) and are nonetheless poorly covered. Conditioning on a description therefore does not explain the gradient away; it sharpens it, from 90 % of described Molecular Key Events down to 50 % of described Population ones, with a test statistic close to twice the unconditional one.

The interpretation is the same but now properly evidenced: the method looks for gene names in prose, so it succeeds where the prose is about molecules and fails where the Key Event is an organ- or organism-level phenotype described in physiological terms. This is the clearest structural limitation of the approach, it is invisible in a single overall percentage, and it is **not** an artefact of uneven curation.


### A5. Coverage by role in the AOP

| role | KEs | described | annotated | % described | % annotated | % annotated of described |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| MIE | 283 | 138 | 123 | 48.8 | 43.5 | 89.1 |
| KE | 1188 | 462 | 357 | 38.9 | 30.1 | 77.3 |
| AO | 238 | 109 | 61 | 45.8 | 25.6 | 56.0 |

**Chi-square test of independence** across roles: chi2(2) = 23.6, p = 7.34e-06 over all Key Events, and chi2(2) = 37.7, p = 6.38e-09 among described Key Events only. Molecular Initiating Events are the best-covered role and Adverse Outcomes the worst — the same gradient as the biological-level result, with the same cause, and again it survives conditioning on the description.


### A6. How concentrated the annotations are

| symbol | hgnc | KEs | % of annotated KEs |
| --- | --- | ---: | ---: |
| ROS1 | 10261 | 32 | 6.6 |
| TNF | 11892 | 29 | 5.9 |
| ESR1 | 3467 | 24 | 4.9 |
| PPARA | 9232 | 24 | 4.9 |
| CD4 | 1678 | 22 | 4.5 |
| COL2A1 | 2200 | 21 | 4.3 |
| IL6 | 6018 | 20 | 4.1 |
| IL1B | 5992 | 19 | 3.9 |
| NFKB1 | 7794 | 18 | 3.7 |
| TGFB1 | 11766 | 17 | 3.5 |
| MAPK1 | 6871 | 17 | 3.5 |
| TP53 | 11998 | 16 | 3.3 |
| MMRN1 | 7178 | 16 | 3.3 |
| IFNA1 | 5417 | 16 | 3.3 |
| PLOD1 | 9081 | 16 | 3.3 |

The ten most frequently detected genes account for **5.2 %** of all KE-gene pairs, while **56.6 %** of genes are detected on exactly one Key Event. The head of this distribution is where dictionary false positives concentrate — short or English-word gene symbols match text that is not about the gene — so the top of this table is the right place to look when judging precision, and it should be inspected before the gene sets are used downstream.


### A6b. Detector disagreement is diagnostic of dictionary false positives

The AOP-Wiki RDF maintainers independently audited the dictionary matcher and documented a set of alias collisions in which an English word or an abbreviation is claimed by a gene symbol (upstream PR #107, open at the time of this release, so its fix is **not** in the pinned data). Those specific symbols can be checked here, and they behave in a way that makes the detector split useful as a quality signal rather than only as a coverage statistic.

| symbol | text probably means | KEs | found by dictionary | found by BERN2 |
| --- | --- | ---: | ---: | ---: |
| ROS1 | reactive oxygen species | 32 | 32 | 0 |
| MMRN1 | extracellular matrix (ECM) | 16 | 15 | 1 |
| TH | thyroid hormone | 14 | 8 | 8 |
| AR | androgen receptor / 'AR' abbreviation | 14 | 6 | 14 |
| TBATA | the adjective 'spatial' | 7 | 7 | 0 |
| AREG | co-asserted with AR | 5 | 5 | 1 |
| FDXR | co-asserted with AR | 5 | 4 | 1 |

**The single most frequently detected gene in the entire dataset is one of them.** ROS1 is annotated onto 32 Key Events, more than any other gene, entirely by the dictionary (32 times) and never once by BERN2. In AOP text 'ROS' is reactive oxygen species. TBATA behaves the same way, and MMRN1 nearly so.

So BERN2 is not merely adding recall — it is **declining** the dictionary's most prominent false positives, which is exactly what a context-aware model should do and what a context-free string match cannot. Two consequences for how this report reads the detector split:

1. The 262 regex-only genes must **not** be presented as 262 genes that would be lost if BERN2 replaced the dictionary. Some are genuine recall the model misses; some are collisions the model correctly rejects. The two cannot be separated without a manually curated benchmark, which does not exist.
2. Regex-only status at high KE frequency is a usable triage signal for curation. The top of the concentration table above is the right place to look for false positives, and it is where the upstream audit found them.


### A7. The Key Event Relationship layer the report does not count

Three of the four free-text fields Method 1 reads belong to KERs, and **703 KERs carry 7752 gene annotations covering 2368 distinct genes** — more genes than the KE layer holds. Section 04 counts only the KE layer.

If KER annotations were propagated to the KEs at either end, a further **111 Key Events** would carry at least one gene, taking coverage from 30.6% to 37.6%. That propagation is not part of Method 1 as defined and is not proposed here — a gene named in the empirical support for a relationship is evidence about the relationship, not necessarily about either event. The figure is reported as an upper bound on what the same text mining already yields but the current KE-level extraction discards.


## B. The organ-specific AOP networks

One frame is reported, for all four organs: the AOP lists Sidra selected and Method 2 was run on, frozen for D3.1 on 2026-08-06. Key Events follow from those AOPs by membership in the pinned release. Two of the 196 liver Key Events differ from the set Method 2 saw, because AOP-Wiki replaced KE 1419 and KE 1420 of AOP 232 with KE 2421 and KE 2423 after her run; the frame follows the release rather than editing a frozen list.


### B1. Coverage and gene counts per organ network

| frame | organ | AOPs | KEs | described KEs | % KEs described | KEs with genes | % KEs with genes | % of described KEs with genes | CI low | CI high | unique genes | KE-gene pairs | median genes / annotated KE | max genes on one KE | genes via regex | genes via BERN2 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| sidra | Brain | 48 | 203 | 118 | 58.1 | 97 | 47.8 | 82.2 | 41.0 | 54.6 | 611 | 986 | 7.0 | 54 | 267 | 490 |
| sidra | Kidney | 16 | 68 | 20 | 29.4 | 14 | 20.6 | 70.0 | 12.7 | 31.6 | 148 | 184 | 11.5 | 29 | 62 | 124 |
| sidra | Liver | 36 | 196 | 91 | 46.4 | 78 | 39.8 | 85.7 | 33.2 | 46.8 | 519 | 811 | 6.0 | 54 | 243 | 415 |
| sidra | Lung | 27 | 134 | 77 | 57.5 | 58 | 43.3 | 75.3 | 35.2 | 51.7 | 626 | 918 | 11.0 | 60 | 290 | 537 |

Both denominators are given. The organ networks are better described than the AOP-Wiki average — 29 % to 58 % of their Key Events carry a description against 39.8 % release-wide — which is unsurprising, since an AOP curated well enough to enter an organ network is curated well enough to have descriptions. So the organ-level coverage figures are **not** comparable with the release-wide 30.6 %, and the report must not present them as though these networks outperform the method. Against their own described populations they run 70–86 %, bracketing the release-wide 76.9 %.


### B2. Regex against BERN2, per organ network

Computed here for the first time — the repository-wide split was known, the per-organ split was not.

| frame | organ | regex only | both | BERN2 only | total genes | % genes BERN2 would have missed | % genes regex would have missed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| sidra | Brain | 121 | 146 | 344 | 611 | 19.8 | 56.3 |
| sidra | Kidney | 24 | 38 | 86 | 148 | 16.2 | 58.1 |
| sidra | Liver | 104 | 139 | 276 | 519 | 20.0 | 53.2 |
| sidra | Lung | 89 | 201 | 336 | 626 | 14.2 | 53.7 |


### B3. Gene overlap between the four organ networks

| region | genes |
| --- | ---: |
| Lung | 307 |
| Brain | 221 |
| Brain & Liver & Lung | 143 |
| Liver | 141 |
| Brain & Liver | 92 |
| Brain & Kidney & Liver & Lung | 68 |
| Brain & Lung | 53 |
| Liver & Lung | 36 |
| Brain & Kidney & Liver | 19 |
| Kidney | 17 |
| Brain & Kidney | 14 |
| Kidney & Liver | 11 |
| Kidney & Lung | 9 |
| Kidney & Liver & Lung | 9 |
| Brain & Kidney & Lung | 1 |

**68 genes are present in all four organ networks.**

Jaccard similarity of the organ gene sets:

|  | Brain | Kidney | Liver | Lung |
| --- | ---: | ---: | ---: | ---: |
| Brain | 1.000 | 0.155 | 0.399 | 0.273 |
| Kidney | 0.155 | 1.000 | 0.191 | 0.127 |
| Liver | 0.399 | 0.191 | 1.000 | 0.288 |
| Lung | 0.273 | 0.127 | 0.288 | 1.000 |

| organ | unique genes | shared with >=1 other | % shared | organ-exclusive |
| --- | ---: | ---: | ---: | ---: |
| Brain | 611 | 390 | 63.8 | 221 |
| Kidney | 148 | 131 | 88.5 | 17 |
| Liver | 519 | 378 | 72.8 | 141 |
| Lung | 626 | 319 | 51.0 | 307 |


### B4. Are the organ networks actually different?

**Coverage** (% of KEs carrying genes), chi-square test of independence across the four networks: chi2(3) = 15.97, p = 0.001. The networks differ in how often their Key Events carry genes.

**Genes per annotated Key Event**, Kruskal-Wallis across the four networks: H(3) = 7.78, p = 0.051. A rank-based test is used because the counts are strongly right-skewed.

Pairwise Mann-Whitney U on genes per annotated KE, Holm-corrected:

| comparison | median A | median B | U | p | p (Holm) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Brain vs Kidney | 7.000 | 11.500 | 469.000 | 0.062 | 0.284 |
| Brain vs Liver | 7.000 | 6.000 | 3750.500 | 0.923 | 1.000 |
| Brain vs Lung | 7.000 | 11.000 | 2229.000 | 0.031 | 0.183 |
| Kidney vs Liver | 11.500 | 6.000 | 713.000 | 0.070 | 0.284 |
| Kidney vs Lung | 11.500 | 11.000 | 426.500 | 0.776 | 1.000 |
| Liver vs Lung | 6.000 | 11.000 | 1829.500 | 0.057 | 0.284 |
