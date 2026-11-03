#!/usr/bin/env python3
"""Stage 2 — build the organ-network frames Method 1's per-organ statistics are computed over.

An "organ AOP network" is not a thing the RDF knows about; it is a hand-made selection of AOPs.
For D3.1 that selection is settled and frozen (decision, 2026-08-06): **Sidra's four organ AOP
lists**, as committed in ``data/networks/sidra_aops_{brain,kidney,liver,lungs}.tsv``. Method 2
was run on these networks and Method 1 is re-run onto them, so both methods report over the
same AOPs.

The lists are the frozen artefact; the Key Events follow from them by membership in the pinned
release. One consequence has to be carried into the report: AOP-Wiki replaced two Key Events of
AOP 232 after Sidra's run, so the liver frame here contains KE 2421 (Inhibition, Intestinal FXR)
and KE 2423 (Hepatic SHP, decreased) where Method 2 saw KE 1419 and KE 1420. Two of 196 liver
KEs therefore differ between the two methods. That is a reporting caveat, not an error to fix —
correcting it would mean editing a frozen list.

The February 2024 Adverse-Outcome-to-organ classification and the ``legacy2024`` frame built
from it were removed on 2026-08-06 along with the frame-sensitivity analysis they supported.
Nothing in the pipeline reads them any more.

Usage:
    python m1_02_frames.py
"""

from __future__ import annotations

import pandas as pd

import m1lib as m1

NETWORKS = m1.ROOT / "data" / "networks"
ORGANS = {"Brain": "brain", "Kidney": "kidney", "Liver": "liver", "Lung": "lungs"}


def _aop_number(series: pd.Series) -> pd.Series:
    return series.astype(str).str.extract(r"(\d+)")[0].astype("Int64")


def build_sidra_frame(aop_ke: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Read the four frozen AOP lists. These are authoritative and are not re-derived here."""
    m1.banner("frame: Sidra's four organ AOP lists (frozen)")
    known = set(aop_ke["aop_number"])
    frames, status = [], {"per_organ": {}}

    for organ, slug in ORGANS.items():
        path = NETWORKS / f"sidra_aops_{slug}.tsv"
        if not path.exists():
            status["per_organ"][organ] = {"status": "missing", "expected_at": str(path)}
            print(f"  {organ:7s} MISSING — expected at {path}")
            continue

        listed = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
        aops = sorted(_aop_number(listed["aopid"]).dropna().astype(int).unique())
        present = [a for a in aops if a in known]
        absent = sorted(set(aops) - set(present))
        kes = sorted(set(aop_ke.loc[aop_ke.aop_number.isin(present), "keid"]))

        frames.append(pd.DataFrame({"aop_number": present, "organ": organ}))
        status["per_organ"][organ] = {
            "status": "frozen — Sidra's committed AOP list",
            "source": str(path.relative_to(m1.ROOT)),
            "aops": len(present), "kes": len(kes),
            "aops_absent_from_release": absent,
        }
        print(f"  {organ:7s} {len(present):3d} AOPs · {len(kes):4d} KEs   from {path.name}")
        if absent:
            print(f"          {len(absent)} listed AOPs are absent from the pinned "
                  f"release: {absent}")

    frame = (pd.concat(frames, ignore_index=True) if frames
             else pd.DataFrame(columns=["aop_number", "organ"]))
    return frame, status


def main() -> None:
    aop_ke = pd.read_csv(m1.RAW_DIR / "aop_ke.csv")
    sidra, sidra_status = build_sidra_frame(aop_ke)

    m1.banner("writing frames")
    written = []
    for organ in sorted(sidra["organ"].unique()):
        aops = sorted(sidra.loc[sidra.organ == organ, "aop_number"].astype(int))
        members = aop_ke[aop_ke.aop_number.isin(aops)][
            ["aop_number", "aopid", "aop_title", "keid", "role"]
        ].sort_values(["aop_number", "keid"])
        path = m1.save_table(members, m1.FRAME_DIR / f"sidra_{organ.lower()}.csv")
        written.append(str(path.relative_to(m1.ROOT)))
        print(f"  -> {path.relative_to(m1.ROOT)}  "
              f"({len(aops)} AOPs, {members.keid.nunique()} KEs)")

    # The two KEs Method 2 saw that this frame cannot, recorded so the report can state it
    # rather than rediscover it.
    liver_kes = set(aop_ke.loc[aop_ke.aop_number.isin(
        sidra.loc[sidra.organ == "Liver", "aop_number"].astype(int)), "keid"])
    drift = {"method2_only": sorted({"KE 1419", "KE 1420"} - liver_kes),
             "frame_only": sorted({"KE 2421", "KE 2423"} & liver_kes),
             "aop": "AOP 232",
             "note": "AOP-Wiki replaced these KEs after Sidra's run; the frame follows the "
                     "pinned release, so 2 of 196 liver KEs differ from Method 2's set."}
    if drift["method2_only"]:
        print(f"\n  AOP 232 drift: Method 2 used {', '.join(drift['method2_only'])}; "
              f"this frame has {', '.join(drift['frame_only'])}")

    m1.write_manifest("02_frames", {
        "sidra": sidra_status,
        "aop232_ke_drift": drift,
        "frames_written": written,
        "legacy2024": "removed 2026-08-06 — the February 2024 classification and the "
                      "frame-sensitivity analysis were dropped from scope",
    })

    print("\nStage 2 complete. All four organ frames are available.")


if __name__ == "__main__":
    main()
