#!/usr/bin/env python3
"""FIG18 — the four organ-specific AOP networks, rendered in Cytoscape.

Reads data/networks/{organ}_nodes.tsv and _edges.tsv (the frozen networks built from
Sidra's AOP selection against the pinned release) and draws each one as Key Events
connected by Key Event Relationships, with MIEs and AOs distinguished from intermediate
Key Events.

Cytoscape must be running with CyREST reachable on http://127.0.0.1:1234. The four
per-organ PNGs land in figures/source/; net_02_assemble.py composes them into the
final panel figure.

Usage:  python analysis/networks/net_01_cytoscape_render.py
"""

from pathlib import Path

import pandas as pd
import py4cytoscape as p4c

ROOT = Path(__file__).resolve().parents[2]
NETDIR = ROOT / "data" / "networks"
OUTDIR = ROOT / "figures" / "source" / "networks"

ORGANS = ["brain", "kidney", "liver", "lungs"]

STYLE = "D31-AOP-network"
COLLECTION = "D3.1 organ AOP networks"

# Role palette. Deliberately not red/green alone — the three roles stay separable in
# greyscale and under the common forms of colour-vision deficiency.
ROLE_COLOURS = {
    "MIE": "#1B7837",  # dark green   — Molecular Initiating Event
    "KE": "#BFD3E6",   # pale blue    — intermediate Key Event
    "AO": "#B2182B",   # dark red     — Adverse Outcome
}
ROLE_SIZES = {"MIE": 45, "KE": 28, "AO": 45}


def load(organ):
    nodes = pd.read_csv(NETDIR / f"{organ}_nodes.tsv", sep="\t", dtype=str)
    edges = pd.read_csv(NETDIR / f"{organ}_edges.tsv", sep="\t", dtype=str)

    # py4cytoscape wants 'id' on nodes and 'source'/'target'/'interaction' on edges.
    nodes = nodes.rename(columns={"name": "label"})
    nodes["role"] = nodes["type"]

    # Drop edges whose endpoints are not in the node table; a KER pointing outside the
    # organ frame would otherwise create a phantom node and inflate the drawn network.
    known = set(nodes["id"])
    before = len(edges)
    edges = edges[edges["source"].isin(known) & edges["target"].isin(known)].copy()
    if len(edges) < before:
        print(f"  {organ}: dropped {before - len(edges)} out-of-frame KERs")

    return nodes, edges


def main():
    p4c.cytoscape_ping()
    print("Cytoscape:", p4c.cytoscape_version_info()["cytoscapeVersion"])
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # A clean style, rebuilt each run so a re-run is deterministic.
    if STYLE in p4c.get_visual_style_names():
        p4c.delete_visual_style(STYLE)
    p4c.create_visual_style(
        STYLE,
        defaults={
            "NODE_SHAPE": "ELLIPSE",
            "NODE_BORDER_WIDTH": 1.0,
            "NODE_BORDER_PAINT": "#444444",
            "NODE_LABEL_FONT_SIZE": 1,
            "NODE_TRANSPARENCY": 255,
            "EDGE_WIDTH": 1.2,
            "EDGE_STROKE_UNSELECTED_PAINT": "#9E9E9E",
            "EDGE_TRANSPARENCY": 160,
            "EDGE_TARGET_ARROW_SHAPE": "DELTA",
            "EDGE_TARGET_ARROW_UNSELECTED_PAINT": "#9E9E9E",
            "NETWORK_BACKGROUND_PAINT": "#FFFFFF",
        },
    )

    for organ in ORGANS:
        print(f"[{organ}]")
        nodes, edges = load(organ)

        suid = p4c.create_network_from_data_frames(
            nodes, edges, title=organ, collection=COLLECTION,
            source_id_list="source", target_id_list="target",
            interaction_type_list="interaction",
        )

        p4c.set_visual_style(STYLE, network=suid)
        p4c.set_node_color_mapping(
            "role", list(ROLE_COLOURS), list(ROLE_COLOURS.values()),
            mapping_type="d", style_name=STYLE, network=suid,
        )
        p4c.set_node_size_mapping(
            "role", list(ROLE_SIZES), list(ROLE_SIZES.values()),
            mapping_type="d", style_name=STYLE, network=suid,
        )
        # cose over force-directed: prefuse strings the disconnected components out in a
        # line beside the giant component, which wastes most of the panel. cose packs them
        # around it, so all four panels come out roughly square and comparable in scale.
        # The layout is stochastic — re-running moves nodes without changing the topology.
        p4c.layout_network("cose", network=suid)
        p4c.fit_content(network=suid)

        out = OUTDIR / f"{organ}.png"
        p4c.export_image(str(out), type="PNG", zoom=250, network=suid,
                         overwrite_file=True)
        print(f"  {len(nodes)} KEs, {len(edges)} KERs -> {out.relative_to(ROOT)}")

    print("\nDone. Run net_02_assemble.py to compose FIG18.")


if __name__ == "__main__":
    main()
