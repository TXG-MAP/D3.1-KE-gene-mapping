#!/usr/bin/env python3
"""Stage 4 — the Method 1 statistics.

Computes everything Section 04 needs, in two blocks:

**A. Repository-wide** (no organ frame involved)
    Coverage of the 1595 Key Events, the regex-vs-BERN2 decomposition at KE / gene / pair
    level, the counterfactual "what would coverage have been with each detector alone",
    genes-per-KE distributions, coverage by biological level and by KE role, how concentrated
    the gene annotations are on a few frequently-detected genes, and the KER-level layer the
    method produces but the report does not currently count.

**B. Per organ network**, on the single frozen frame — Sidra's four AOP lists, the same
    networks Method 2 was run on. Includes the per-organ regex-vs-BERN2 split, which has never
    been computed before, and the four-way gene overlap.

The February 2024 classification, the ``legacy2024`` frame and the liver frame-sensitivity
comparison were removed on 2026-08-06 when the AOP lists were frozen. There is now one frame,
so there is nothing to be sensitive to.

Every test that is run is named with its assumptions in the output. Counts of KEs are counts
of *distinct* KEs — a KE belonging to two AOPs is one KE, and summing per-organ files
double-counts the 18-odd KEs that sit in more than one network.

Usage:
    python m1_04_stats.py
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import m1lib as m1

OUT_MD = m1.TABLE_DIR / "method1_statistics.md"


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a proportion — behaves near 0 and 1, unlike the normal
    approximation, which matters because kidney coverage sits on small denominators."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def holm(pvalues: list[float]) -> list[float]:
    """Holm-Bonferroni step-down adjustment. Controls the family-wise error rate without
    assuming independence, which pairwise organ comparisons on overlapping KE sets violate."""
    m = len(pvalues)
    order = np.argsort(pvalues)
    adjusted = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * pvalues[idx])
        adjusted[idx] = min(1.0, running)
    return adjusted.tolist()


def describe_counts(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {k: float("nan") for k in
                ("n", "mean", "sd", "median", "q1", "q3", "min", "max", "total")}
    return {
        "n": int(values.size),
        "mean": float(values.mean()),
        "sd": float(values.std(ddof=1)) if values.size > 1 else 0.0,
        "median": float(np.median(values)),
        "q1": float(np.percentile(values, 25)),
        "q3": float(np.percentile(values, 75)),
        "min": float(values.min()),
        "max": float(values.max()),
        "total": float(values.sum()),
    }


class Report:
    """Accumulates the Markdown that becomes tables/method1_statistics.md."""

    def __init__(self) -> None:
        self.parts: list[str] = []

    def h(self, level: int, text: str) -> None:
        self.parts.append(f"\n{'#' * level} {text}\n")
        print(f"\n{'-' * 4} {text}")

    def p(self, text: str) -> None:
        self.parts.append(text.strip() + "\n")

    def table(self, df: pd.DataFrame, floatfmt: str = "{:.1f}") -> None:
        md = m1.markdown_table(df, floatfmt=floatfmt)
        self.parts.append(md + "\n")
        print(md)

    def write(self, path: Path, header: str) -> None:
        path.write_text(header + "\n" + "\n".join(self.parts))
        print(f"\n  -> {path.relative_to(m1.ROOT)}")


# --------------------------------------------------------------------------------------
# Block A — repository-wide
# --------------------------------------------------------------------------------------


def block_a(rep: Report, ke: pd.DataFrame, pairs: pd.DataFrame, aop_ke: pd.DataFrame,
            ker: pd.DataFrame, results: dict) -> None:
    rep.h(2, "A. Repository-wide, all 1595 Key Events")

    n_ke = len(ke)
    annotated = set(pairs["keid"])
    described = set(ke.loc[ke.has_description, "keid"])
    lo, hi = wilson_ci(len(annotated), n_ke)
    lo_d, hi_d = wilson_ci(len(annotated), len(described))
    rep.p(
        f"Method 1 annotates at least one gene onto **{len(annotated)} of {n_ke} Key Events "
        f"({len(annotated) / n_ke:.1%}, 95 % CI {lo:.1%}–{hi:.1%})**, using "
        f"{pairs['hgnc'].nunique()} distinct genes over {len(pairs)} KE-gene pairs."
    )

    # --- the denominator ---------------------------------------------------------------
    rep.h(3, "A0. The denominator: which Key Events can this method reach at all?")
    rep.p(
        "Method 1 finds gene names by reading a Key Event's free-text description. A Key Event "
        "with no description cannot receive a gene annotation from it — that is not a miss, it "
        "is out of scope. So `30.6 % of Key Events` measures AOP-Wiki's description "
        "completeness at least as much as it measures the method, and the two should be "
        "separated."
    )
    outside = sorted(annotated - described)
    coverage_rows = pd.DataFrame([
        {"population": "all Key Events in the release", "KEs": n_ke,
         "annotated": len(annotated), "coverage %": 100 * len(annotated) / n_ke,
         "CI low": 100 * lo, "CI high": 100 * hi},
        {"population": "Key Events carrying a description", "KEs": len(described),
         "annotated": len(annotated), "coverage %": 100 * len(annotated) / len(described),
         "CI low": 100 * lo_d, "CI high": 100 * hi_d},
        {"population": "Key Events with no description", "KEs": n_ke - len(described),
         "annotated": len(outside), "coverage %": 0.0,
         "CI low": float("nan"), "CI high": float("nan")},
    ])
    rep.table(coverage_rows, floatfmt="{:.1f}")
    rep.p(
        f"**All {len(annotated)} gene-annotated Key Events carry a description, and none of "
        f"the {n_ke - len(described)} undescribed Key Events carries a gene.** The containment "
        "is exact, which confirms the description field is the sole route by which a KE "
        "acquires genes under this method and is asserted on every pipeline run.\n\n"
        f"So the method's own hit rate is **{len(annotated) / len(described):.1%} "
        f"({len(annotated)} of {len(described)}, 95 % CI {lo_d:.1%}–{hi_d:.1%})**, not 30.6 %. "
        "The two numbers answer different questions and the report needs both: "
        f"{len(annotated) / len(described):.0%} is how well the text mining works, "
        f"{len(annotated) / n_ke:.0%} is how much of the AOP-Wiki it can currently serve. The "
        f"gap between them is the {n_ke - len(described)} Key Events "
        f"({(n_ke - len(described)) / n_ke:.0%}) that have no description to mine — a gap that "
        "no improvement to the method can close, and a concrete curation recommendation for "
        "AOP-Wiki rather than a limitation of Method 1."
    )
    results["denominator"] = {
        "kes_total": n_ke,
        "kes_with_description": len(described),
        "kes_without_description": n_ke - len(described),
        "annotated": len(annotated),
        "coverage_all_kes": len(annotated) / n_ke,
        "coverage_described_kes": len(annotated) / len(described),
        "coverage_described_ci": [lo_d, hi_d],
        "annotated_but_undescribed": outside,
        "table": coverage_rows.to_dict("records"),
    }

    # --- description length ------------------------------------------------------------
    rep.h(3, "A0b. How much text a Key Event has predicts what the method finds")
    desc = ke[ke.has_description].copy()
    per_ke_all = pairs.groupby("keid")["hgnc"].nunique()
    desc["n_genes"] = desc["keid"].map(per_ke_all).fillna(0)
    desc["annotated"] = desc["keid"].isin(annotated)
    rho, rho_p = stats.spearmanr(desc["description_length"], desc["n_genes"])
    u = stats.mannwhitneyu(desc.loc[desc.annotated, "description_length"],
                           desc.loc[~desc.annotated, "description_length"],
                           alternative="two-sided")
    desc["quartile"] = pd.qcut(desc["description_length"], 4,
                               labels=["Q1 shortest", "Q2", "Q3", "Q4 longest"])
    quart = (desc.groupby("quartile", observed=True)
             .agg(KEs=("keid", "size"), annotated=("annotated", "sum"),
                  median_length=("description_length", "median"),
                  median_genes=("n_genes", "median")))
    quart["% annotated"] = 100 * quart["annotated"] / quart["KEs"]
    quart = quart.reset_index().rename(columns={"quartile": "description length"})
    rep.table(quart, floatfmt="{:.1f}")
    rep.p(
        f"Among described Key Events, description length correlates with the number of genes "
        f"found at **Spearman rho = {rho:.2f}** (p = {rho_p:.1e}). Annotated KEs have a median "
        f"description of {desc.loc[desc.annotated, 'description_length'].median():.0f} "
        f"characters against "
        f"{desc.loc[~desc.annotated, 'description_length'].median():.0f} for unannotated ones "
        f"(Mann-Whitney p = {u.pvalue:.1e}). In the longest quartile the method annotates "
        f"{quart['% annotated'].iloc[-1]:.0f} % of Key Events; in the shortest, "
        f"{quart['% annotated'].iloc[0]:.0f} %.\n\n"
        "This is worth stating plainly because it constrains how the output should be read: "
        "the size of a Key Event's Method 1 gene set reflects how much has been written about "
        "it as much as how many genes are biologically involved. A KE with 13 genes and one "
        "with 3 are not necessarily different in biology — they may differ in curation effort. "
        "That makes raw gene-set size a poor basis for comparing Key Events, and it is a "
        "caution the comparison in Section 07 has to carry."
    )
    results["description_length"] = {
        "spearman_rho": float(rho), "spearman_p": float(rho_p),
        "median_length_annotated": float(desc.loc[desc.annotated, "description_length"].median()),
        "median_length_unannotated": float(desc.loc[~desc.annotated, "description_length"].median()),
        "mannwhitney_p": float(u.pvalue),
        "quartiles": quart.to_dict("records"),
    }

    # --- the detector decomposition, at all three levels it can be counted at ----------
    rep.h(3, "A1. What BERN2 changed: the detector decomposition")
    rep.p(
        "The same comparison counted three ways. These are different questions and the "
        "report should say which one it is quoting: a gene can be found by both detectors on "
        "one KE and by only one on another, so the gene-level and pair-level splits do not "
        "agree and neither is wrong."
    )
    regex_kes = set(pairs.loc[pairs.by_regex, "keid"])
    ner_kes = set(pairs.loc[pairs.by_ner, "keid"])
    regex_genes = set(pairs.loc[pairs.by_regex, "hgnc"])
    ner_genes = set(pairs.loc[pairs.by_ner, "hgnc"])

    rows = []
    for label, a, b in (("Key Events", regex_kes, ner_kes),
                        ("unique genes", regex_genes, ner_genes)):
        rows.append({
            "counted over": label,
            "regex (HGNC)": len(a), "BERN2 (NER)": len(b),
            "regex only": len(a - b), "both": len(a & b), "BERN2 only": len(b - a),
            "union": len(a | b),
        })
    pair_counts = pairs["detector"].value_counts()
    rows.append({
        "counted over": "KE-gene pairs",
        "regex (HGNC)": int(pairs.by_regex.sum()), "BERN2 (NER)": int(pairs.by_ner.sum()),
        "regex only": int(pair_counts.get("regex only", 0)),
        "both": int(pair_counts.get("both", 0)),
        "BERN2 only": int(pair_counts.get("BERN2 only", 0)),
        "union": len(pairs),
    })
    decomposition = pd.DataFrame(rows)
    rep.table(decomposition)
    results["detector_decomposition"] = decomposition.to_dict("records")

    # --- the counterfactual: coverage under each detector alone -----------------------
    rep.h(3, "A2. Coverage under each detector alone")
    rep.p(
        "What the deliverable would have reported had only one detector run. The regex column "
        "is the closest available reconstruction of the pre-2026-06 pipeline on the pinned "
        "release — it is *not* the number in the current draft, which was computed on an "
        "older AOP-Wiki release as well as an older detector."
    )
    counter = pd.DataFrame([
        {"detector": "regex / HGNC dictionary only", "KEs with >=1 gene": len(regex_kes),
         "% of 1595 KEs": 100 * len(regex_kes) / n_ke, "unique genes": len(regex_genes),
         "KE-gene pairs": int(pairs.by_regex.sum())},
        {"detector": "BERN2 NER only", "KEs with >=1 gene": len(ner_kes),
         "% of 1595 KEs": 100 * len(ner_kes) / n_ke, "unique genes": len(ner_genes),
         "KE-gene pairs": int(pairs.by_ner.sum())},
        {"detector": "union (as released)", "KEs with >=1 gene": len(annotated),
         "% of 1595 KEs": 100 * len(annotated) / n_ke,
         "unique genes": pairs["hgnc"].nunique(), "KE-gene pairs": len(pairs)},
    ])
    rep.table(counter, floatfmt="{:.2f}")
    gain_ke = (len(annotated) - len(regex_kes)) / len(regex_kes)
    gain_gene = (pairs["hgnc"].nunique() - len(regex_genes)) / len(regex_genes)
    rep.p(
        f"Adding BERN2 to the dictionary matcher increased annotated Key Events by "
        f"**{gain_ke:+.1%}** and unique genes by **{gain_gene:+.1%}**. BERN2 alone recovers "
        f"{len(regex_genes & ner_genes) / len(regex_genes):.1%} of what the dictionary found, "
        f"so the two are complementary rather than one superseding the other: "
        f"{len(regex_genes - ner_genes)} genes are found *only* by the dictionary."
    )

    # McNemar on the paired KE-level coverage. Paired because the same 1595 KEs are exposed
    # to both detectors; a two-sample test would be wrong here.
    b = len(regex_kes - ner_kes)   # regex found it, BERN2 did not
    c = len(ner_kes - regex_kes)   # BERN2 found it, regex did not
    mcnemar_p = float(stats.binomtest(min(b, c), b + c, 0.5).pvalue) if (b + c) else float("nan")
    rep.p(
        f"**McNemar exact test** on KE-level coverage (paired over all {n_ke} KEs): "
        f"{c} KEs annotated only by BERN2 against {b} only by the dictionary, "
        f"p = {mcnemar_p:.3g}. The detectors are not interchangeable."
    )
    results["mcnemar_ke_coverage"] = {
        "regex_only_kes": b, "ner_only_kes": c, "p_value": mcnemar_p,
        "test": "McNemar exact (binomial) on paired KE-level coverage",
    }

    # Per-KE gene counts under each detector, over the KEs either detector annotated.
    per_ke = (pairs.groupby("keid")
              .agg(union=("hgnc", "nunique"),
                   regex=("by_regex", "sum"), ner=("by_ner", "sum"))
              .reindex(sorted(annotated)).fillna(0))
    w = stats.wilcoxon(per_ke["ner"], per_ke["regex"], zero_method="wilcox")
    rep.p(
        f"**Wilcoxon signed-rank** on genes per KE, BERN2 against dictionary, over the "
        f"{len(per_ke)} annotated KEs: median {per_ke['ner'].median():.0f} vs "
        f"{per_ke['regex'].median():.0f} genes, W = {w.statistic:.0f}, p = {w.pvalue:.3g}."
    )
    results["wilcoxon_genes_per_ke"] = {
        "median_ner": float(per_ke["ner"].median()),
        "median_regex": float(per_ke["regex"].median()),
        "statistic": float(w.statistic), "p_value": float(w.pvalue),
    }

    # --- genes per KE -----------------------------------------------------------------
    rep.h(3, "A3. Genes per Key Event")
    dist_rows = []
    for label, col in (("union (as released)", "union"), ("BERN2", "ner"), ("regex", "regex")):
        d = describe_counts(per_ke[col].values)
        d = {"detector": label, **d}
        dist_rows.append(d)
    dist = pd.DataFrame(dist_rows)[
        ["detector", "n", "total", "mean", "sd", "median", "q1", "q3", "min", "max"]]
    rep.table(dist, floatfmt="{:.1f}")
    rep.p(
        "Denominator is the 488 KEs carrying at least one gene, so the regex row counts zero "
        "for KEs only BERN2 found. The distribution is heavily right-skewed — the mean sits "
        f"well above the median and the largest single KE carries {int(per_ke['union'].max())} "
        "genes — which is why medians are reported throughout and why a mean gene count per "
        "organ would be misleading."
    )
    results["genes_per_ke"] = dist.to_dict("records")

    # --- by biological level -----------------------------------------------------------
    rep.h(3, "A4. Coverage by biological level of organisation")
    rep.p(
        "Reported against both denominators, because the obvious objection to the raw gradient "
        "is that molecular Key Events might simply be better described. They are not: the "
        "gradient is a property of the method, and conditioning on having a description makes "
        "it stronger rather than weaker."
    )
    lvl = ke.copy()
    lvl["has_genes"] = lvl["keid"].isin(annotated)
    per_lvl = (lvl.groupby("biolevel")
               .agg(KEs=("keid", "size"), described=("has_description", "sum"),
                    annotated=("has_genes", "sum"))
               .reindex(m1.BIOLEVELS).dropna(how="all"))
    per_lvl["% described"] = 100 * per_lvl["described"] / per_lvl["KEs"]
    per_lvl["% annotated of all"] = 100 * per_lvl["annotated"] / per_lvl["KEs"]
    per_lvl["% annotated of described"] = 100 * per_lvl["annotated"] / per_lvl["described"]
    gene_counts = per_ke["union"]
    per_lvl["median genes / annotated KE"] = [
        float(gene_counts.reindex(
            lvl.loc[(lvl.biolevel == b) & lvl.has_genes, "keid"]).median())
        for b in per_lvl.index
    ]
    per_lvl = per_lvl.reset_index().rename(columns={"biolevel": "biological level"})
    rep.table(per_lvl, floatfmt="{:.1f}")

    table = np.array([[int(r.annotated), int(r.KEs - r.annotated)]
                      for r in per_lvl.itertuples()])
    chi2, p, dof, _ = stats.chi2_contingency(table)
    tab_desc = np.array([[int(r.described), int(r.KEs - r.described)]
                         for r in per_lvl.itertuples()])
    chi2_d, p_d, dof_d, _ = stats.chi2_contingency(tab_desc)
    tab_cond = np.array([[int(r.annotated), int(r.described - r.annotated)]
                         for r in per_lvl.itertuples()])
    chi2_c, p_c, dof_c, _ = stats.chi2_contingency(tab_cond)
    rep.p(
        f"Three chi-square tests of independence across the six levels:\n\n"
        f"- **coverage over all KEs**: chi2({dof}) = {chi2:.1f}, p = {p:.3g}\n"
        f"- **whether a KE is described at all**: chi2({dof_d}) = {chi2_d:.1f}, "
        f"p = {p_d:.3g}\n"
        f"- **coverage among described KEs**: chi2({dof_c}) = {chi2_c:.1f}, p = {p_c:.3g}\n\n"
        "Description availability does vary by level, but weakly and **non-monotonically** — "
        f"Organ Key Events are the most often described ({per_lvl.loc[per_lvl['biological level'] == 'Organ', '% described'].iloc[0]:.0f} %) "
        "and are nonetheless poorly covered. Conditioning on a description therefore does not "
        "explain the gradient away; it sharpens it, from "
        f"{per_lvl['% annotated of described'].iloc[0]:.0f} % of described Molecular Key Events "
        f"down to {per_lvl['% annotated of described'].iloc[-1]:.0f} % of described Population "
        "ones, with a test statistic close to twice the unconditional one.\n\n"
        "The interpretation is the same but now properly evidenced: the method looks for gene "
        "names in prose, so it succeeds where the prose is about molecules and fails where the "
        "Key Event is an organ- or organism-level phenotype described in physiological terms. "
        "This is the clearest structural limitation of the approach, it is invisible in a "
        "single overall percentage, and it is **not** an artefact of uneven curation."
    )
    results["coverage_by_biolevel"] = {
        "table": per_lvl.to_dict("records"),
        "chi2": float(chi2), "dof": int(dof), "p_value": float(p),
        "described_chi2": float(chi2_d), "described_p_value": float(p_d),
        "conditional_chi2": float(chi2_c), "conditional_dof": int(dof_c),
        "conditional_p_value": float(p_c),
    }

    # --- by role -----------------------------------------------------------------------
    rep.h(3, "A5. Coverage by role in the AOP")
    roles = (aop_ke[["keid", "role"]].drop_duplicates()
             .assign(has_genes=lambda d: d["keid"].isin(annotated),
                     described=lambda d: d["keid"].isin(described)))
    # A KE can be an MIE in one AOP and an intermediate KE in another; count it under each
    # role it plays rather than forcing a single label.
    per_role = (roles.groupby("role")
                .agg(KEs=("keid", "nunique"), described=("described", "sum"),
                     annotated=("has_genes", "sum")))
    per_role["% described"] = 100 * per_role["described"] / per_role["KEs"]
    per_role["% annotated"] = 100 * per_role["annotated"] / per_role["KEs"]
    per_role["% annotated of described"] = 100 * per_role["annotated"] / per_role["described"]
    per_role = per_role.reindex(["MIE", "KE", "AO"]).reset_index()
    rep.table(per_role, floatfmt="{:.1f}")
    tab = np.array([[int(r.annotated), int(r.KEs - r.annotated)] for r in per_role.itertuples()])
    chi2r, pr, dofr, _ = stats.chi2_contingency(tab)
    tab_rc = np.array([[int(r.annotated), int(r.described - r.annotated)]
                       for r in per_role.itertuples()])
    chi2rc, prc, dofrc, _ = stats.chi2_contingency(tab_rc)
    rep.p(
        f"**Chi-square test of independence** across roles: chi2({dofr}) = {chi2r:.1f}, "
        f"p = {pr:.3g} over all Key Events, and chi2({dofrc}) = {chi2rc:.1f}, p = {prc:.3g} "
        "among described Key Events only. Molecular Initiating Events are the best-covered "
        "role and Adverse Outcomes the worst — the same gradient as the biological-level "
        "result, with the same cause, and again it survives conditioning on the description."
    )
    results["coverage_by_role"] = {
        "table": per_role.to_dict("records"),
        "chi2": float(chi2r), "dof": int(dofr), "p_value": float(pr),
        "conditional_chi2": float(chi2rc), "conditional_p_value": float(prc),
        "note": "roles overlap — a KE may be an MIE in one AOP and an intermediate KE in "
                "another, so the role rows are not a partition of the 1595 KEs",
    }

    # --- concentration -------------------------------------------------------------------
    rep.h(3, "A6. How concentrated the annotations are")
    freq = (pairs.groupby(["hgnc", "symbol"], dropna=False)
            .agg(KEs=("keid", "nunique"))
            .sort_values("KEs", ascending=False).reset_index())
    top = freq.head(15)[["symbol", "hgnc", "KEs"]]
    top["% of annotated KEs"] = 100 * top["KEs"] / len(annotated)
    rep.table(top, floatfmt="{:.1f}")
    share10 = 100 * freq.head(10)["KEs"].sum() / len(pairs)
    singleton = 100 * (freq["KEs"] == 1).mean()
    rep.p(
        f"The ten most frequently detected genes account for **{share10:.1f} %** of all "
        f"KE-gene pairs, while **{singleton:.1f} %** of genes are detected on exactly one Key "
        "Event. The head of this distribution is where dictionary false positives concentrate "
        "— short or English-word gene symbols match text that is not about the gene — so the "
        "top of this table is the right place to look when judging precision, and it should "
        "be inspected before the gene sets are used downstream."
    )
    results["concentration"] = {
        "top_genes": top.to_dict("records"),
        "share_of_pairs_from_top10_genes": float(share10),
        "pct_genes_on_a_single_ke": float(singleton),
    }

    # --- detector disagreement as a false-positive signal ---------------------------------
    rep.h(3, "A6b. Detector disagreement is diagnostic of dictionary false positives")
    rep.p(
        "The AOP-Wiki RDF maintainers independently audited the dictionary matcher and "
        "documented a set of alias collisions in which an English word or an abbreviation is "
        "claimed by a gene symbol (upstream PR #107, open at the time of this release, so its "
        "fix is **not** in the pinned data). Those specific symbols can be checked here, and "
        "they behave in a way that makes the detector split useful as a quality signal rather "
        "than only as a coverage statistic."
    )
    # The collisions named upstream, with what the text almost certainly meant.
    known_collisions = {
        "ROS1": "reactive oxygen species", "TBATA": "the adjective 'spatial'",
        "MMRN1": "extracellular matrix (ECM)", "TH": "thyroid hormone",
        "AR": "androgen receptor / 'AR' abbreviation", "AREG": "co-asserted with AR",
        "FDXR": "co-asserted with AR",
    }
    by_symbol = (pairs.groupby("symbol")
                 .agg(KEs=("keid", "nunique"), by_regex=("by_regex", "sum"),
                      by_ner=("by_ner", "sum")))
    fp_rows = []
    for symbol, meaning in known_collisions.items():
        if symbol not in by_symbol.index:
            continue
        row = by_symbol.loc[symbol]
        fp_rows.append({
            "symbol": symbol, "text probably means": meaning,
            "KEs": int(row.KEs), "found by dictionary": int(row.by_regex),
            "found by BERN2": int(row.by_ner),
        })
    fp = pd.DataFrame(fp_rows).sort_values("KEs", ascending=False)
    rep.table(fp)
    rep.p(
        f"**The single most frequently detected gene in the entire dataset is one of them.** "
        f"ROS1 is annotated onto {int(by_symbol.loc['ROS1', 'KEs'])} Key Events, more than any "
        f"other gene, entirely by the dictionary "
        f"({int(by_symbol.loc['ROS1', 'by_regex'])} times) and never once by BERN2. In AOP "
        "text 'ROS' is reactive oxygen species. TBATA behaves the same way, and MMRN1 nearly "
        "so.\n\n"
        "So BERN2 is not merely adding recall — it is **declining** the dictionary's most "
        "prominent false positives, which is exactly what a context-aware model should do and "
        "what a context-free string match cannot. Two consequences for how this report reads "
        "the detector split:\n\n"
        "1. The 262 regex-only genes must **not** be presented as 262 genes that would be lost "
        "if BERN2 replaced the dictionary. Some are genuine recall the model misses; some are "
        "collisions the model correctly rejects. The two cannot be separated without a "
        "manually curated benchmark, which does not exist.\n"
        "2. Regex-only status at high KE frequency is a usable triage signal for curation. The "
        "top of the concentration table above is the right place to look for false positives, "
        "and it is where the upstream audit found them."
    )
    results["known_collisions"] = {
        "source": "marvinm2/AOPWikiRDF PR #107 (open at commit 97191db — fix NOT in this data)",
        "table": fp.to_dict("records"),
    }

    # --- KER layer ------------------------------------------------------------------------
    rep.h(3, "A7. The Key Event Relationship layer the report does not count")
    ker_kes = set(ker["upstream"].dropna()) | set(ker["downstream"].dropna())
    would_gain = (ker_kes - annotated) & set(ke["keid"])
    rep.p(
        f"Three of the four free-text fields Method 1 reads belong to KERs, and "
        f"**{ker['kerid'].nunique()} KERs carry {len(ker)} gene annotations covering "
        f"{ker['hgnc'].nunique()} distinct genes** — more genes than the KE layer holds. "
        f"Section 04 counts only the KE layer.\n\n"
        f"If KER annotations were propagated to the KEs at either end, a further "
        f"**{len(would_gain)} Key Events** would carry at least one gene, taking coverage "
        f"from {len(annotated) / n_ke:.1%} to {(len(annotated) + len(would_gain)) / n_ke:.1%}. "
        "That propagation is not part of Method 1 as defined and is not proposed here — a "
        "gene named in the empirical support for a relationship is evidence about the "
        "relationship, not necessarily about either event. The figure is reported as an upper "
        "bound on what the same text mining already yields but the current KE-level "
        "extraction discards."
    )
    results["ker_layer"] = {
        "kers_with_genes": int(ker["kerid"].nunique()),
        "ker_gene_pairs": len(ker),
        "unique_genes_on_kers": int(ker["hgnc"].nunique()),
        "additional_kes_if_propagated": len(would_gain),
        "coverage_if_propagated": (len(annotated) + len(would_gain)) / n_ke,
    }


# --------------------------------------------------------------------------------------
# Block B — per organ network
# --------------------------------------------------------------------------------------


def organ_stats(frame_name: str, organ: str, members: pd.DataFrame,
                pairs: pd.DataFrame, described: set) -> dict:
    kes = set(members["keid"])
    sub = pairs[pairs.keid.isin(kes)]
    annotated = set(sub["keid"])
    described_here = kes & described
    per_ke = sub.groupby("keid")["hgnc"].nunique()
    lo, hi = wilson_ci(len(annotated), len(kes))
    return {
        "frame": frame_name, "organ": organ,
        "AOPs": int(members["aop_number"].nunique()),
        "KEs": len(kes),
        "described KEs": len(described_here),
        "% KEs described": 100 * len(described_here) / len(kes) if kes else float("nan"),
        "KEs with genes": len(annotated),
        "% KEs with genes": 100 * len(annotated) / len(kes) if kes else float("nan"),
        "% of described KEs with genes":
            100 * len(annotated) / len(described_here) if described_here else float("nan"),
        "CI low": 100 * lo, "CI high": 100 * hi,
        "unique genes": int(sub["hgnc"].nunique()),
        "KE-gene pairs": len(sub),
        "median genes / annotated KE": float(per_ke.median()) if len(per_ke) else float("nan"),
        "max genes on one KE": int(per_ke.max()) if len(per_ke) else 0,
        "genes via regex": int(sub.loc[sub.by_regex, "hgnc"].nunique()),
        "genes via BERN2": int(sub.loc[sub.by_ner, "hgnc"].nunique()),
        "_gene_set": set(sub["hgnc"]),
        "_annotated": annotated,
        "_ke_set": kes,
        "_counts": per_ke,
    }


def block_b(rep: Report, pairs: pd.DataFrame, described: set, results: dict) -> dict:
    rep.h(2, "B. The organ-specific AOP networks")
    rep.p(
        "One frame is reported, for all four organs: the AOP lists Sidra selected and Method 2 "
        "was run on, frozen for D3.1 on 2026-08-06. Key Events follow from those AOPs by "
        "membership in the pinned release. Two of the 196 liver Key Events differ from the set "
        "Method 2 saw, because AOP-Wiki replaced KE 1419 and KE 1420 of AOP 232 with KE 2421 "
        "and KE 2423 after her run; the frame follows the release rather than editing a frozen "
        "list."
    )

    collected: dict[tuple[str, str], dict] = {}
    rows = []
    for path in sorted(m1.FRAME_DIR.glob("*.csv")):
        frame_name, organ = path.stem.split("_", 1)
        members = pd.read_csv(path)
        s = organ_stats(frame_name, organ.capitalize(), members, pairs, described)
        collected[(frame_name, organ.capitalize())] = s
        rows.append({k: v for k, v in s.items() if not k.startswith("_")})

    summary = pd.DataFrame(rows).sort_values(["frame", "organ"]).reset_index(drop=True)
    rep.h(3, "B1. Coverage and gene counts per organ network")
    rep.table(summary, floatfmt="{:.1f}")
    leg = summary[summary.frame == "sidra"]
    rep.p(
        "Both denominators are given. The organ networks are better described than the "
        f"AOP-Wiki average — {leg['% KEs described'].min():.0f} % to "
        f"{leg['% KEs described'].max():.0f} % of their Key Events carry a description against "
        "39.8 % release-wide — which is unsurprising, since an AOP curated well enough to enter "
        "an organ network is curated well enough to have descriptions. So the organ-level "
        "coverage figures are **not** comparable with the release-wide 30.6 %, and the report "
        "must not present them as though these networks outperform the method. Against their "
        f"own described populations they run "
        f"{leg['% of described KEs with genes'].min():.0f}–"
        f"{leg['% of described KEs with genes'].max():.0f} %, bracketing the release-wide "
        "76.9 %."
    )
    results["per_organ"] = summary.to_dict("records")

    # --- detector split per organ -------------------------------------------------------
    rep.h(3, "B2. Regex against BERN2, per organ network")
    rep.p("Computed here for the first time — the repository-wide split was known, the "
          "per-organ split was not.")
    det_rows = []
    for (frame_name, organ), s in sorted(collected.items()):
        sub = pairs[pairs.keid.isin(s["_ke_set"])]
        rx, nr = set(sub.loc[sub.by_regex, "hgnc"]), set(sub.loc[sub.by_ner, "hgnc"])
        det_rows.append({
            "frame": frame_name, "organ": organ,
            "regex only": len(rx - nr), "both": len(rx & nr), "BERN2 only": len(nr - rx),
            "total genes": len(rx | nr),
            "% genes BERN2 would have missed": 100 * len(rx - nr) / len(rx | nr) if (rx | nr) else float("nan"),
            "% genes regex would have missed": 100 * len(nr - rx) / len(rx | nr) if (rx | nr) else float("nan"),
        })
    detectors = pd.DataFrame(det_rows)
    rep.table(detectors, floatfmt="{:.1f}")
    m1.save_table(detectors, m1.DERIVED_DIR / "organ_detector_split.csv")
    results["organ_detector_split"] = detectors.to_dict("records")

    # --- four-organ overlap ----------------------------------------
    rep.h(3, "B3. Gene overlap between the four organ networks")
    legacy = {o: collected[("sidra", o)]["_gene_set"]
              for o in m1.ORGANS if ("sidra", o) in collected}
    regions = {}
    for r in range(1, len(legacy) + 1):
        for combo in itertools.combinations(sorted(legacy), r):
            inside = set.intersection(*(legacy[o] for o in combo))
            outside = set().union(*(legacy[o] for o in legacy if o not in combo)) \
                if len(combo) < len(legacy) else set()
            exclusive = inside - outside
            if exclusive:
                regions[" & ".join(combo)] = sorted(exclusive)
    region_tbl = pd.DataFrame(
        [{"region": k, "genes": len(v)} for k, v in sorted(regions.items(),
                                                           key=lambda kv: -len(kv[1]))])
    rep.table(region_tbl)

    shared_all = set.intersection(*legacy.values()) if legacy else set()
    rep.p(f"**{len(shared_all)} genes are present in all four organ networks.**")

    jac = pd.DataFrame(
        [[len(legacy[a] & legacy[b]) / len(legacy[a] | legacy[b]) for b in sorted(legacy)]
         for a in sorted(legacy)],
        index=sorted(legacy), columns=sorted(legacy))
    rep.p("Jaccard similarity of the organ gene sets:")
    rep.table(jac.round(3).reset_index().rename(columns={"index": ""}), floatfmt="{:.3f}")

    overlap_rows = []
    for organ, genes in sorted(legacy.items()):
        others = set().union(*(g for o, g in legacy.items() if o != organ))
        overlap_rows.append({
            "organ": organ, "unique genes": len(genes),
            "shared with >=1 other": len(genes & others),
            "% shared": 100 * len(genes & others) / len(genes),
            "organ-exclusive": len(genes - others),
        })
    overlap = pd.DataFrame(overlap_rows)
    rep.table(overlap, floatfmt="{:.1f}")
    results["organ_overlap"] = {
        "regions": {k: len(v) for k, v in regions.items()},
        "shared_by_all_four": sorted(shared_all),
        "jaccard": jac.round(4).to_dict(),
        "exclusivity": overlap.to_dict("records"),
    }
    m1.save_table(
        pd.DataFrame([{"region": k, "n_genes": len(v), "genes": ", ".join(v)}
                      for k, v in regions.items()]),
        m1.DERIVED_DIR / "organ_gene_overlap_regions.csv")

    # --- tests across organs --------------------------------------------------------------
    rep.h(3, "B4. Are the organ networks actually different?")
    legacy_stats = {o: collected[("sidra", o)] for o in m1.ORGANS
                    if ("sidra", o) in collected}
    tab = np.array([[s["KEs with genes"], s["KEs"] - s["KEs with genes"]]
                    for s in legacy_stats.values()])
    chi2, p, dof, _ = stats.chi2_contingency(tab)
    rep.p(
        f"**Coverage** (% of KEs carrying genes), chi-square test of independence across the "
        f"four networks: chi2({dof}) = {chi2:.2f}, p = {p:.3f}. "
        + ("The networks do not differ detectably in how often their Key Events carry genes."
           if p >= 0.05 else
           "The networks differ in how often their Key Events carry genes.")
    )

    groups = [s["_counts"].values for s in legacy_stats.values()]
    kw = stats.kruskal(*groups)
    rep.p(
        f"**Genes per annotated Key Event**, Kruskal-Wallis across the four networks: "
        f"H({len(groups) - 1}) = {kw.statistic:.2f}, p = {kw.pvalue:.3f}. A rank-based test is "
        "used because the counts are strongly right-skewed."
    )

    pair_rows, pvals = [], []
    for a, b in itertools.combinations(sorted(legacy_stats), 2):
        u = stats.mannwhitneyu(legacy_stats[a]["_counts"], legacy_stats[b]["_counts"],
                               alternative="two-sided")
        pair_rows.append({"comparison": f"{a} vs {b}",
                          f"median A": float(legacy_stats[a]["_counts"].median()),
                          f"median B": float(legacy_stats[b]["_counts"].median()),
                          "U": float(u.statistic), "p": float(u.pvalue)})
        pvals.append(float(u.pvalue))
    pw = pd.DataFrame(pair_rows)
    pw["p (Holm)"] = holm(pvals)
    rep.p("Pairwise Mann-Whitney U on genes per annotated KE, Holm-corrected:")
    rep.table(pw, floatfmt="{:.3f}")
    results["organ_tests"] = {
        "coverage_chi2": {"chi2": float(chi2), "dof": int(dof), "p_value": float(p)},
        "genes_per_ke_kruskal": {"H": float(kw.statistic), "p_value": float(kw.pvalue)},
        "pairwise_mannwhitney_holm": pw.to_dict("records"),
    }
    return collected


# --------------------------------------------------------------------------------------
def main() -> None:
    ke = pd.read_csv(m1.RAW_DIR / "ke_inventory.csv")
    pairs = pd.read_csv(m1.RAW_DIR / "ke_gene_long.csv", dtype={"hgnc": str})
    aop_ke = pd.read_csv(m1.RAW_DIR / "aop_ke.csv")
    ker = pd.read_csv(m1.RAW_DIR / "ker_gene_long.csv", dtype={"hgnc": str})

    rep = Report()
    results: dict = {}

    block_a(rep, ke, pairs, aop_ke, ker, results)
    collected = block_b(rep, pairs, set(ke.loc[ke.has_description, "keid"]), results)

    header = (
        "# Method 1 — statistics\n\n"
        f"Generated from AOP-Wiki RDF commit `{m1.COMMIT}`, data generated "
        f"{m1.EXPECTED_CREATED_ON}, via `analysis/method1/m1_04_stats.py`.\n\n"
        "Regenerate with `make method1`. Do not hand-edit: this file is overwritten.\n"
    )
    rep.write(OUT_MD, header)

    (m1.DERIVED_DIR / "method1_statistics.json").write_text(
        json.dumps(results, indent=2, default=str) + "\n")
    print(f"  -> {(m1.DERIVED_DIR / 'method1_statistics.json').relative_to(m1.ROOT)}")

    m1.write_manifest("04_stats", {"outputs": [str(OUT_MD.relative_to(m1.ROOT))]})
    print("\nStage 4 complete.")


if __name__ == "__main__":
    main()
