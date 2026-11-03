"""Shared configuration, SPARQL plumbing and plot style for the Method 1 analysis.

Everything the numbered ``m1_0*.py`` scripts have in common lives here so that no script
hard-codes an endpoint, a release date or a colour. The endpoint and the pin are read from
``pinned-release.yaml`` at the repository root, which is the single machine-readable record of
what D3.1 is pinned to.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]

PINNED_RELEASE = ROOT / "pinned-release.yaml"
QUERY_DIR = ROOT / "queries"

RAW_DIR = ROOT / "data" / "method1" / "raw"
FRAME_DIR = ROOT / "data" / "method1" / "frames"
DERIVED_DIR = ROOT / "data" / "method1" / "derived"
TABLE_DIR = ROOT / "tables"
FIG_FINAL = ROOT / "figures" / "final"
FIG_SOURCE = ROOT / "figures" / "source"
FIG_CLEAN = ROOT / "figures" / "clean"

for _d in (RAW_DIR, FRAME_DIR, DERIVED_DIR, TABLE_DIR, FIG_FINAL, FIG_SOURCE):
    _d.mkdir(parents=True, exist_ok=True)

# Caption-free build. With D31_FIG_CLEAN=1 the provenance strip under each plot is dropped
# and the figures are written to figures/clean/ instead of figures/final/, leaving the
# captioned versions untouched. The report carries that text as the figure caption, so the
# image should not repeat it; figures/final/ stays the self-contained standalone version.
CLEAN_FIGURES = os.environ.get("D31_FIG_CLEAN", "") == "1"

# The 2025-02 regex-era outputs sit directly in data/method1/ and are historical evidence for
# what the current draft says. Nothing here may overwrite them.
LEGACY_OUTPUTS = {
    ROOT / "data" / "method1" / f"{o}_table.csv" for o in ("Brain", "Kidney", "Liver", "Lung")
} | {ROOT / "data" / "method1" / "Intersections_4.txt"}


# --------------------------------------------------------------------------------------
# The pin
# --------------------------------------------------------------------------------------


def load_pin() -> dict:
    """Return the ``aopwiki_rdf`` block of pinned-release.yaml (the canonical endpoint)."""
    with open(PINNED_RELEASE) as fh:
        return yaml.safe_load(fh)["aopwiki_rdf"]


PIN = load_pin()
ENDPOINT = PIN["endpoint"]
EXPECTED_CREATED_ON = str(PIN["pav_created_on"])
EXPECTED_VOID = {k: int(v) for k, v in PIN["void_triples"].items()}
COMMIT = PIN["commit"]


# --------------------------------------------------------------------------------------
# SPARQL
# --------------------------------------------------------------------------------------


def _sparql_client():
    from SPARQLWrapper import JSON, SPARQLWrapper

    client = SPARQLWrapper(ENDPOINT)
    client.setReturnFormat(JSON)
    client.setTimeout(600)
    return client


def run_query(query: str, retries: int = 4, pause: float = 3.0) -> pd.DataFrame:
    """Execute a SPARQL SELECT and return the bindings as a DataFrame of plain strings.

    Unbound OPTIONAL variables become NA rather than being dropped, so a column that is
    absent from every row still exists in the result.
    """
    client = _sparql_client()
    client.setQuery(query)
    last = None
    for attempt in range(1, retries + 1):
        try:
            result = client.query().convert()
            break
        except Exception as exc:  # network flakiness against a shared endpoint
            last = exc
            if attempt == retries:
                raise RuntimeError(f"SPARQL query failed after {retries} attempts: {exc}") from exc
            print(f"  ! query attempt {attempt} failed ({exc}); retrying in {pause:.0f}s",
                  file=sys.stderr)
            time.sleep(pause * attempt)
    else:  # pragma: no cover - defensive
        raise RuntimeError(f"SPARQL query failed: {last}")

    cols = result["head"]["vars"]
    rows = [{c: binding.get(c, {}).get("value") for c in cols}
            for binding in result["results"]["bindings"]]
    return pd.DataFrame(rows, columns=cols)


def read_query(name: str) -> str:
    """Load a ``.rq`` file from queries/. Comment lines are kept — they are valid SPARQL."""
    path = QUERY_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"query file not found: {path}")
    return path.read_text()


def verify_pin(allow_drift: bool = False) -> dict:
    """Check the live endpoint still serves the pinned release.

    The endpoint regenerates weekly and will drift off the pinned commit. ``pav:createdOn``
    catches a different generation date; the three ``void:triples`` counts distinguish two
    runs made on the same date. Returns the observed values.
    """
    created = run_query(
        "PREFIX pav: <http://purl.org/pav/> SELECT DISTINCT ?c WHERE { ?s pav:createdOn ?c }"
    )
    void = run_query(
        "PREFIX void: <http://rdfs.org/ns/void#> SELECT ?s ?t WHERE { ?s void:triples ?t }"
    )
    observed_created = sorted(created["c"].dropna().unique().tolist())
    observed_void = {
        row["s"].rsplit("/", 1)[-1]: int(row["t"]) for _, row in void.iterrows()
    }
    # The main graph is registered without the .ttl suffix in pav but with it in void.
    observed_void.pop("AOPWikiRDF", None)

    problems = []
    if observed_created != [EXPECTED_CREATED_ON]:
        problems.append(
            f"pav:createdOn is {observed_created}, expected ['{EXPECTED_CREATED_ON}']")
    for key, expected in EXPECTED_VOID.items():
        got = observed_void.get(key)
        if got != expected:
            problems.append(f"void:triples[{key}] is {got}, expected {expected}")

    status = {
        "endpoint": ENDPOINT,
        "pinned_commit": COMMIT,
        "expected_created_on": EXPECTED_CREATED_ON,
        "observed_created_on": observed_created,
        "expected_void_triples": EXPECTED_VOID,
        "observed_void_triples": observed_void,
        "matches_pin": not problems,
        "problems": problems,
    }

    if problems:
        message = (
            "ENDPOINT HAS DRIFTED OFF THE PIN:\n  "
            + "\n  ".join(problems)
            + "\n\nThe endpoint is rolling and regenerates weekly. Either re-pin "
              "pinned-release.yaml to the release now being served, or load the TTLs by SHA "
              f"from the {COMMIT} data_dumps URLs. Re-run with --allow-drift to proceed "
              "anyway (results will NOT be reproducible against the stated pin)."
        )
        if not allow_drift:
            raise SystemExit(message)
        print("WARNING: " + message, file=sys.stderr)
    else:
        print(f"  pin verified: {EXPECTED_CREATED_ON}, commit {COMMIT}, void counts match")
    return status


# --------------------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------------------


def git_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], text=True
        ).strip()
    except Exception:
        return "unknown"


def write_manifest(name: str, payload: dict) -> Path:
    """Write a provenance manifest next to the outputs of one pipeline stage."""
    payload = dict(payload)
    payload.update({
        "stage": name,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "aopwiki_rdf_commit": COMMIT,
        "aopwiki_rdf_data_generated": EXPECTED_CREATED_ON,
        "endpoint": ENDPOINT,
        "d31_repo_revision": git_revision(),
        "python": sys.version.split()[0],
    })
    path = RAW_DIR.parent / f"manifest_{name}.json"
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    return path


def save_table(df: pd.DataFrame, path: Path, index: bool = False) -> Path:
    """Write a CSV, refusing to clobber the 2025-02 regex-era artefacts."""
    path = Path(path)
    if path.resolve() in {p.resolve() for p in LEGACY_OUTPUTS if p.exists()}:
        raise RuntimeError(
            f"refusing to overwrite the archived regex-era artefact {path}. "
            "Those files are the evidence for what the current draft says."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index)
    return path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------
# Domain constants
# --------------------------------------------------------------------------------------

ORGANS = ["Brain", "Kidney", "Liver", "Lung"]

#: Biological level of organisation, ordered from molecular upwards rather than alphabetically.
BIOLEVELS = ["Molecular", "Cellular", "Tissue", "Organ", "Individual", "Population"]

#: Detector provenance of a single (KE, gene) pair.
DETECTORS = ["regex only", "both", "BERN2 only"]

#: Organ-network frames. `legacy2024` is the February 2024 AO->organ classification the
#: current draft was computed on; `sidra` is the authoritative selection strategy (Q9), of
#: which only liver has been written down. See m1_02_frames.py.
FRAMES = ["legacy2024", "sidra"]


# --------------------------------------------------------------------------------------
# Plot style
# --------------------------------------------------------------------------------------

#: Categorical hues for the four organ networks, in fixed order — never cycled, never
#: reassigned when a subset is plotted. Checked for colour-vision-deficiency separation (worst CVD deltaE 9.2, worst normal-vision deltaE 16.3, light mode):
#: all checks pass. The aqua sits just under 3:1 against the surface, so every figure that
#: uses it also carries direct labels or a table view, which is the documented relief.
ORGAN_COLOURS = {
    "Brain": "#2a78d6",   # blue
    "Kidney": "#eb6834",  # orange
    "Liver": "#1baf7a",   # aqua
    "Lung": "#4a3aa7",    # violet
}

#: Detector provenance, likewise validated under `--pairs all` (worst CVD deltaE 9.2).
DETECTOR_COLOURS = {
    "regex only": "#eb6834",
    "both": "#1baf7a",
    "BERN2 only": "#4a3aa7",
}

#: Single-hue sequential ramp, light->dark, validated as ordinal (monotone lightness, visible
#: step gaps, light end clears 2:1 against the surface).
SEQUENTIAL = ["#7fb0e0", "#4f8fd4", "#2a6cc4", "#17457f"]

SURFACE = "#ffffff"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8880"
GRID = "#e3e2dd"


def apply_style() -> None:
    """Recessive grid, thin marks, text in ink tokens rather than series colour."""
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlepad": 8,
        "axes.labelsize": 9,
        "axes.labelcolor": INK_SECONDARY,
        "axes.edgecolor": GRID,
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": INK_SECONDARY,
        "ytick.color": INK_SECONDARY,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "text.color": INK,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "lines.linewidth": 2.0,
        "lines.markersize": 5,
    })


def despine(ax, keep=("left", "bottom")) -> None:
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def source_note(fig, y: float, text: str, **kwargs):
    """The muted provenance strip under a plot — release, frame, and the test statistics.

    Suppressed under D31_FIG_CLEAN=1. The y offsets are negative and the figures save with
    bbox_inches="tight", so dropping the note also drops the canvas it was extending into.
    """
    if CLEAN_FIGURES:
        return None
    kwargs.setdefault("fontsize", 8)
    kwargs.setdefault("color", INK_MUTED)
    kwargs.setdefault("ha", "left")
    return fig.text(0.008, y, text, **kwargs)


def save_figure(fig, stem: str, final: bool = True) -> Path:
    """Write a figure as PNG (report) and PDF (vector, for the archive)."""
    target_dir = FIG_CLEAN if CLEAN_FIGURES else (FIG_FINAL if final else FIG_SOURCE)
    target_dir.mkdir(parents=True, exist_ok=True)
    png = target_dir / f"{stem}.png"
    fig.savefig(png)
    fig.savefig(target_dir / f"{stem}.pdf")
    print(f"  -> {png.relative_to(ROOT)}")
    return png


def fmt_pct(x: float, digits: int = 1) -> str:
    return f"{x:.{digits}f} %"


def markdown_table(df: pd.DataFrame, floatfmt: str = "{:.1f}") -> str:
    """Render a DataFrame as a GitHub-flavoured pipe table, numeric columns right-aligned."""
    def cell(v):
        if isinstance(v, float):
            return "" if pd.isna(v) else floatfmt.format(v)
        return "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v)

    cols = list(df.columns)
    numeric = [pd.api.types.is_numeric_dtype(df[c]) for c in cols]
    head = "| " + " | ".join(cols) + " |"
    rule = "| " + " | ".join("---:" if n else "---" for n in numeric) + " |"
    body = ["| " + " | ".join(cell(v) for v in row) + " |"
            for row in df.itertuples(index=False, name=None)]
    return "\n".join([head, rule, *body])


def banner(text: str) -> None:
    print(f"\n=== {text} " + "=" * max(0, 72 - len(text)))
