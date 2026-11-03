"""Shared paths and provenance helpers for Method 2's gene-level layer.

Method 2's curation (KE -> GO Biological Process) is published on Zenodo
(concept DOI 10.5281/zenodo.21264952). What lives here is only the step the draft
described but never implemented: GO term -> gene set, harmonised to Ensembl so the three
methods can be compared.

``m1lib`` owns the AOP-Wiki pin, the plot style and the checksum helper, so those are imported
rather than copied. Its ``write_manifest`` is NOT reused: it writes into ``data/method1/``
regardless of caller, so a Method 2 stage calling it silently files its provenance under
Method 1. The local one below writes next to Method 2's own outputs.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method1"))
import m1lib as m1  # noqa: E402
from m1lib import COMMIT, EXPECTED_CREATED_ON, banner, git_revision, sha256  # noqa: E402,F401

ROOT = m1.ROOT
M2 = ROOT / "data" / "method2"
CACHE = M2 / "cache"
DERIVED = M2 / "derived"

for _d in (M2, CACHE, DERIVED):
    _d.mkdir(parents=True, exist_ok=True)

CURATED = M2 / "curated_ke_go_mapping.csv"

FOUR = ["brain", "kidney", "liver", "lung"]
NODE_FILE = {"brain": "brain_nodes.tsv", "kidney": "kidney_nodes.tsv",
             "liver": "liver_nodes.tsv", "lung": "lungs_nodes.tsv"}


def write_manifest(name: str, payload: dict) -> Path:
    payload = dict(payload)
    payload.update({
        "stage": name,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "aopwiki_rdf_commit": COMMIT,
        "aopwiki_rdf_data_generated": EXPECTED_CREATED_ON,
        "d31_repo_revision": git_revision(),
        "python": sys.version.split()[0],
    })
    path = M2 / f"manifest_{name}.json"
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print(f"  -> {path.relative_to(ROOT)}")
    return path


def save_table(df, path: Path, index: bool = False, sep: str = "\t") -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index, sep=sep)
    print(f"  -> {path.relative_to(ROOT)}  ({len(df)} rows)")
    return path


def load_networks() -> dict[str, set[str]]:
    """The four frozen organ networks, as sets of KE ids."""
    import pandas as pd
    return {
        organ: set(pd.read_csv(ROOT / "data" / "networks" / NODE_FILE[organ],
                               sep="\t")["id"].str.strip())
        for organ in FOUR
    }
