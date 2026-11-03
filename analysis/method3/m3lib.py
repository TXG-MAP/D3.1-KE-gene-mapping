"""Shared configuration, HTTP plumbing and scoring for the Method 3 (KE->WikiPathways) run.

Method 3 reads ranked pathway suggestions from the molAOP Builder and adjudicates them
offline. Nothing here writes to the service: every call is a GET, and the two write verbs the
builder exposes (``POST /submit``, ``POST /admin/.../approve``) are deliberately absent.

The AOP-Wiki pin is read from ``pinned-release.yaml`` through ``m1lib``, so Methods 1 and 3
cannot drift onto different releases. The builder's own corpus has no pin upstream; stage 2
creates one by snapshotting it and recording the export revision.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# Method 1 owns the pin, the provenance helpers and the plot style. Import rather than copy so
# a re-pin lands on both methods at once.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "method1"))
import m1lib  # noqa: E402
from m1lib import COMMIT, EXPECTED_CREATED_ON, banner, git_revision, sha256  # noqa: E402,F401

ROOT = m1lib.ROOT

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------

M3 = ROOT / "data" / "method3"
FRAME_DIR = M3 / "frames"
CORPUS_DIR = M3 / "corpus"          # tracked: the offline copy of what was queried
CACHE_DIR = M3 / "cache" / "suggest"  # gitignored: raw per-query ranker responses
DERIVED_DIR = M3 / "derived"
ADJ_DIR = M3 / "adjudication"
SHEET_DIR = M3 / "sheets"           # TSV twin of every workbook sheet
INPUTS_DIR = Path(__file__).resolve().parent / "inputs"

for _d in (FRAME_DIR, CORPUS_DIR, CACHE_DIR, DERIVED_DIR, ADJ_DIR, SHEET_DIR, INPUTS_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------------------
# The organ networks
# --------------------------------------------------------------------------------------

#: The three networks this run covers. Liver was curated in July 2026 and is out of scope.
ORGANS = ["brain", "kidney", "lung"]

#: Sidra's node/edge files spell lung "lungs"; the report and Method 1 say "lung".
NODE_FILE = {"brain": "brain_nodes.tsv", "kidney": "kidney_nodes.tsv", "lung": "lungs_nodes.tsv"}

#: Expected KE counts, asserted in stage 1. These come from the frozen networks (Q9, settled
#: 2026-08-06) and are the definition of scope — a mismatch means the frames moved.
EXPECTED_KES = {"brain": 203, "kidney": 68, "lung": 134}
EXPECTED_UNION = 375

# --------------------------------------------------------------------------------------
# The builder
# --------------------------------------------------------------------------------------

import os  # noqa: E402

BASE = os.environ.get("MOLAOP_BASE", "https://molaop-builder.vhp4safety.nl").rstrip("/")

#: /suggest_pathways, /search_pathways, /get_ke_options and /get_pathway_options share one
#: bucket: 500 requests/hour (src/services/rate_limiter.py:193 at builder HEAD af422df).
SPARQL_BUCKET_PER_HOUR = 500

#: Seconds between requests. 500/h is one every 7.2 s; this run stays inside a window per organ
#: rather than pacing to the cap, so 8 s leaves headroom for the snapshot calls.
PACE = 8.0

USER_AGENT = "d31-method3/1.0 (EFSA TXG-MAP; read-only)"

#: Incremented by every request that actually leaves the machine. Cache hits do not count.
REQUEST_COUNT = 0


class RateLimited(SystemExit):
    """Raised on a 429. Never retried — the limiter records rejected calls too, so retrying
    deepens the hole. Resume from the cache after the window rolls over."""


def get(path: str, params: dict | None = None, timeout: int = 90,
        pace: float | None = None, raw: bool = False):
    """One GET against the builder. Returns parsed JSON, or (bytes, headers) when ``raw``.

    5xx is transient — the single CPU-bound replica drops one through Traefik now and then —
    and is retried with backoff. 429 is fatal by design.
    """
    global REQUEST_COUNT
    url = f"{BASE}{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                REQUEST_COUNT += 1
                body = resp.read()
                # Header names are case-insensitive on the wire but dict() preserves the
                # server's casing, so `X-Export-Revision` misses a lookup for the lowercase
                # name and the provenance field silently records "".
                headers = {k.lower(): v for k, v in resp.headers.items()}
            time.sleep(PACE if pace is None else pace)
            return (body, headers) if raw else json.loads(body)
        except urllib.error.HTTPError as exc:
            REQUEST_COUNT += 1
            if exc.code == 429:
                raise RateLimited(
                    f"\n429 on {url}\n"
                    f"The limiter counts rejected calls, so retrying makes it worse. Wait for "
                    f"the window to roll over ({SPARQL_BUCKET_PER_HOUR}/h) and re-run — the "
                    f"cache resumes where this stopped."
                ) from exc
            if exc.code in (500, 502, 503, 504) and attempt < 3:
                time.sleep(10 * (attempt + 1))
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if attempt < 3:
                time.sleep(10 * (attempt + 1))
                continue
            raise
    raise RuntimeError(f"unreachable: {url}")


# --------------------------------------------------------------------------------------
# The builder's confidence formula
# --------------------------------------------------------------------------------------
#
# Authority: molAOP-builder/scoring_config.yaml block `ke_pathway_assessment:` and
# src/core/assessment_scoring.py. Reproduced here so the tier can be predicted offline, and
# validated against every live approved mapping by m3_00_validate_predictor.py before any
# adjudication runs. step1 is NOT scored; it only supplies the connection type.

EVIDENCE = {"known": 3, "likely": 2, "possible": 1, "uncertain": 0}
SPECIFICITY = {"specific": 2, "includes": 1, "loose": 0}
COVERAGE = {"complete": 1.5, "keysteps": 1.0, "minor": 0.5}

#: The bonus test is a case-insensitive SUBSTRING match on the biological level, per
#: assessment_scoring.py's `resolve_biolevel`. "Cellular" and "cellular level" both qualify.
BONUS_LEVELS = ("molecular", "cellular", "tissue")
BONUS = 1.0

HIGH, MEDIUM = 5.0, 2.5

CONNECTION = {
    "causative": "causative",
    "responsive": "responsive",
    "bidirectional": "other",
    "unclear": "undefined",
}

STEP_OPTIONS = {
    "step1": tuple(CONNECTION),
    "step2": tuple(EVIDENCE),
    "step3": tuple(SPECIFICITY),
    "step4": tuple(COVERAGE),
}


def biolevel_bonus(biolevel: str | None) -> float:
    level = (biolevel or "").lower()
    return BONUS if any(q in level for q in BONUS_LEVELS) else 0.0


def predict(step2: str, step3: str, step4: str, biolevel: str | None) -> tuple[float, str]:
    """Return (score, tier). An unanswered question contributes 0 rather than aborting."""
    score = (EVIDENCE.get(step2, 0)
             + SPECIFICITY.get(step3, 0)
             + COVERAGE.get(step4, 0)
             + biolevel_bonus(biolevel))
    tier = "high" if score >= HIGH else "medium" if score >= MEDIUM else "low"
    return round(score, 2), tier


def would_submit(step3: str, step4: str, tier: str) -> tuple[bool, str]:
    """Apply the pre-registered not-submitting rules. Nothing here submits anything — this
    only records what a curator following the rubric would lodge."""
    if step3 == "loose" and step4 == "minor":
        return False, "loose + minor is a rejection, not a Low submission"
    if tier == "low":
        return False, "predicted low; an honest empty slot beats a weak mapping"
    return True, ""


# --------------------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------------------


def write_manifest(name: str, payload: dict) -> Path:
    """Write a provenance manifest for one Method 3 stage, next to the outputs."""
    payload = dict(payload)
    payload.update({
        "stage": name,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "aopwiki_rdf_commit": COMMIT,
        "aopwiki_rdf_data_generated": EXPECTED_CREATED_ON,
        "molaop_builder_base": BASE,
        "d31_repo_revision": git_revision(),
        "python": sys.version.split()[0],
        "requests_issued": REQUEST_COUNT,
    })
    path = M3 / f"manifest_{name}.json"
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print(f"  -> {path.relative_to(ROOT)}")
    return path


def save_table(df, path: Path, index: bool = False, sep: str = "\t") -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index, sep=sep)
    print(f"  -> {path.relative_to(ROOT)}  ({len(df)} rows)")
    return path


def load_json(path: Path):
    with open(path) as fh:
        return json.load(fh)


def dump_json(obj, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=str) + "\n")
    print(f"  -> {path.relative_to(ROOT)}")
    return path
