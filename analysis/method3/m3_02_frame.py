#!/usr/bin/env python3
"""Stage 2 — build the Key Event frame for the brain, kidney and lung networks.

One row per KE, carrying everything the rest of the pipeline keys on. The frame is the
definition of scope, so it asserts its own counts: Sidra's frozen networks hold 203 / 68 / 134
KEs and 375 in the union, and a drift in any of those means the networks moved under the run
rather than that a KE was dropped.

Three title columns, deliberately:

  ke_title          the pinned release's dc:title — the authoritative one for the report
  ke_title_node     the string in Sidra's node table, kept so a divergence is visible
  ke_title_builder  what /get_ke_options serves, which is what the GUI sends to the ranker

The raw query arm must use `ke_title_builder`, not the release title, or the run stops being
reproducible by a reviewer clicking through the interface — which is the single most important
caveat Method 3 has to state.

Offline: reads only what stage 1 captured.
"""
import csv
import json
from collections import defaultdict

import pandas as pd

import m3lib as L


def read_nodes() -> pd.DataFrame:
    """KE identity and MIE/KE/AO role, per organ, from the frozen networks."""
    frames = []
    for organ in L.ORGANS:
        path = L.ROOT / "data" / "networks" / L.NODE_FILE[organ]
        df = pd.read_csv(path, sep="\t")
        if len(df) != L.EXPECTED_KES[organ]:
            raise SystemExit(
                f"{path.name} holds {len(df)} KEs, expected {L.EXPECTED_KES[organ]}. The organ "
                "networks were frozen on 2026-08-06 and are not to be "
                "regenerated — reconcile before continuing.")
        df["organ"] = organ
        frames.append(df.rename(columns={"id": "ke_id", "name": "ke_title_node", "type": "role"}))
    return pd.concat(frames, ignore_index=True)


def read_aop_context() -> dict[str, list[str]]:
    """AOP membership per KE, from Sidra's authoritative AOP->KE frames."""
    ctx = defaultdict(set)
    for organ in L.ORGANS:
        path = L.ROOT / "data" / "method1" / "frames" / f"sidra_{organ}.csv"
        for row in csv.DictReader(open(path)):
            ctx[row["keid"]].add(row["aopid"])
    return {k: sorted(v, key=lambda a: int(a.split()[1])) for k, v in ctx.items()}


def main():
    L.banner("stage 2: build the KE frame for brain, kidney and lung")

    nodes = read_nodes()
    aop_ctx = read_aop_context()

    release = pd.read_csv(L.CORPUS_DIR / "ke_descriptions.csv").fillna(
        {"description": "", "biolevel": "", "title": ""})
    rel = release.set_index("keid")

    builder = {k["KElabel"]: k for k in json.load(open(L.CORPUS_DIR / "ke_options.json"))}

    approved = defaultdict(list)
    for m in json.load(open(L.CORPUS_DIR / "live_mappings.json"))["data"]:
        approved[m["ke_id"]].append(m["pathway_id"])

    # Collapse to one row per KE, with organ membership as a column. 50 KEs sit in more than
    # one of the four networks, so a per-organ row layout would double-count them in every
    # downstream total.
    organs_of = nodes.groupby("ke_id")["organ"].apply(lambda s: ";".join(sorted(set(s))))
    roles_of = nodes.groupby("ke_id")["role"].apply(lambda s: ";".join(sorted(set(s))))
    node_title = nodes.groupby("ke_id")["ke_title_node"].first()

    rows = []
    for ke_id in sorted(organs_of.index, key=lambda x: int(x.split()[1])):
        r = rel.loc[ke_id] if ke_id in rel.index else None
        b = builder.get(ke_id)
        desc = (r["description"] if r is not None else "") or ""
        biolevel = (r["biolevel"] if r is not None else "") or ""
        rows.append({
            "ke_id": ke_id,
            "ke_title": (r["title"] if r is not None else "") or node_title[ke_id],
            "ke_title_node": node_title[ke_id],
            "ke_title_builder": (b or {}).get("KEtitle", ""),
            "organs": organs_of[ke_id],
            "n_organs": len(organs_of[ke_id].split(";")),
            "role": roles_of[ke_id],
            "biolevel": biolevel,
            "biolevel_bonus": L.biolevel_bonus(biolevel),
            "aop_ids": ";".join(aop_ctx.get(ke_id, [])),
            "n_aops": len(aop_ctx.get(ke_id, [])),
            "has_description": int(bool(desc)),
            "description_length": len(desc),
            "in_builder": int(b is not None),
            "n_approved_wp": len(approved.get(ke_id, [])),
            "approved_wps": ";".join(sorted(approved.get(ke_id, []))),
        })

    frame = pd.DataFrame(rows)

    if len(frame) != L.EXPECTED_UNION:
        raise SystemExit(f"frame holds {len(frame)} KEs, expected {L.EXPECTED_UNION}")

    # A KE whose builder title differs from the release title will be ranked on a string the
    # report does not quote. Not fatal, but it has to be visible.
    drift = frame[(frame["in_builder"] == 1)
                  & (frame["ke_title_builder"].str.strip() != frame["ke_title"].str.strip())]

    print(f"\n  {len(frame)} distinct KEs across brain, kidney and lung")
    for organ in L.ORGANS:
        sub = frame[frame["organs"].str.contains(organ)]
        print(f"    {organ:7s} {len(sub):3d} KEs · {sub['has_description'].sum():3d} described "
              f"· {int(sub['biolevel_bonus'].sum() / L.BONUS):3d} earn the biolevel bonus "
              f"· {(sub['n_approved_wp'] > 0).sum():3d} already mapped")
    print(f"    shared  {(frame['n_organs'] > 1).sum():3d} KEs sit in more than one of the three")
    print(f"\n  missing from the builder: {(frame['in_builder'] == 0).sum()}")
    print(f"  builder title differs from the release title: {len(drift)}")
    print(f"  biolevels: {frame['biolevel'].value_counts().to_dict()}")

    L.save_table(frame, L.FRAME_DIR / "ke_frame.tsv")
    if len(drift):
        L.save_table(drift[["ke_id", "ke_title", "ke_title_builder", "organs"]],
                     L.FRAME_DIR / "ke_title_drift.tsv")

    L.write_manifest("02_frame", {
        "kes": len(frame),
        "per_organ": {o: int(frame["organs"].str.contains(o).sum()) for o in L.ORGANS},
        "shared_across_networks": int((frame["n_organs"] > 1).sum()),
        "described": int(frame["has_description"].sum()),
        "earning_biolevel_bonus": int((frame["biolevel_bonus"] > 0).sum()),
        "already_mapped": int((frame["n_approved_wp"] > 0).sum()),
        "absent_from_builder": int((frame["in_builder"] == 0).sum()),
        "builder_title_drift": len(drift),
    })


if __name__ == "__main__":
    main()
