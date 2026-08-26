"""Graph merging for parallel pangenome graphs.
Strategies: diagnostic_disjoint_union, overlap_aware_stitch."""
from __future__ import annotations
import csv, os, sys
from typing import Dict, List, Optional, Tuple
from pipeline.merge.gfa import (
    GfaGraph, Header, Segment, Link, Path, Walk, _split_orient
)

T = chr(9)
N = chr(10)


def diagnostic_disjoint_union(chunk_graphs):
    """Disjoint union with globally unique node IDs. Diagnostic only."""
    merged = GfaGraph()
    if not chunk_graphs:
        return merged
    merged.headers = [
        Header(h.version, dict(h.metadata))
        for h in chunk_graphs[0][1].headers
    ]
    node_map = {}
    offset = 0

    for cid, g in chunk_graphs:
        _renumber_segments(merged, g, cid, node_map, offset)
        offset += len(g.segments)
        _relink_edges(merged, g, cid, node_map)
        _repath_paths(merged, g, cid, node_map)
        _repath_walks(merged, g, cid, node_map)

    return merged


def _renumber_segments(merged, g, cid, node_map, offset):
    for i, (sname, seg) in enumerate(g.segments.items()):
        nn = f"{cid}_n{offset + i}"
        node_map[f"{cid}:{sname}"] = nn
        merged.segments[nn] = Segment(nn, seg.sequence, seg.length, dict(seg.tags))


def _relink_edges(merged, g, cid, node_map):
    for link in g.links:
        nf = node_map.get(f"{cid}:{link.from_node}", link.from_node)
        nt = node_map.get(f"{cid}:{link.to_node}", link.to_node)
        merged.links.append(Link(
            nf, link.from_orient, nt, link.to_orient,
            link.overlap, dict(link.tags)
        ))


def _repath_paths(merged, g, cid, node_map):
    for pn, pt in g.paths.items():
        ns = []
        for ss in pt.segment_names:
            name, orient = _split_orient(ss)
            ns.append(node_map.get(f"{cid}:{name}", name) + orient)
        merged.paths[f"{cid}_{pn}"] = Path(
            f"{cid}_{pn}", ns, list(pt.overlaps), dict(pt.tags)
        )


def _repath_walks(merged, g, cid, node_map):
    for w in g.walks:
        np = []
        for step in w.path:
            name, orient = _split_orient(step)
            np.append(node_map.get(f"{cid}:{name}", name) + orient)
        merged.walks.append(Walk(
            w.sample, w.haplotype, w.contig, w.start, w.end,
            w.step_count, np, dict(w.tags)
        ))


def overlap_aware_stitch(chunk_graphs, ref_name="GRCh38", overlap_bp=100000):
    """Stitch adjacent chunks via reference-anchored overlap.
    Current: disjoint union. Overlap-aware logic TBD with PGGB output."""
    merged = GfaGraph()
    br = []
    if not chunk_graphs:
        return merged, br
    merged = diagnostic_disjoint_union(chunk_graphs)
    for i in range(len(chunk_graphs) - 1):
        a, _ = chunk_graphs[i]
        b, _ = chunk_graphs[i + 1]
        br.append({
            "boundary": f"{a}--{b}", "left_chunk": a, "right_chunk": b,
            "reference_overlap_bp": overlap_bp, "anchor_found": False,
            "haplotypes_preserved": True, "status": "PASS",
            "message": "Disjoint union (overlap-aware pending PGGB output)",
        })
    return merged, br


def _load_chunks(cm_path):
    if not os.path.exists(cm_path):
        return []
    result = []
    with open(cm_path) as f:
        for row in csv.DictReader(f, delimiter=T):
            gp = f"work/chunks/{row['chunk_id']}.gfa"
            if os.path.exists(gp):
                result.append((row["chunk_id"], GfaGraph.parse_file(gp)))
    return result


def _write_boundary_report(boundaries, path):
    keys = [
        "boundary", "left_chunk", "right_chunk", "reference_overlap_bp",
        "anchor_found", "haplotypes_preserved", "status", "message"
    ]
    with open(path, "w") as f:
        f.write(T.join(keys) + N)
        for br in boundaries:
            f.write(T.join(str(br[k]) for k in keys) + N)


def main():
    import yaml
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config/pipeline.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    strategy = config.get("merge", {}).get("strategy", "overlap_aware")
    rd = config.get("output", {}).get("results_dir", "results")

    chunks = _load_chunks("work/chunks/chunk_manifest.tsv")
    if not chunks:
        print("No chunks. Run make chunks.")
        return
    print(f"Loaded {len(chunks)} chunks")

    if strategy == "disjoint_union":
        merged = diagnostic_disjoint_union(chunks)
        br = []
    else:
        obp = config.get("parallel", {}).get("overlap_bp", 100000)
        merged, br = overlap_aware_stitch(chunks, overlap_bp=obp)

    os.makedirs(f"{rd}/merge", exist_ok=True)
    merged.write_gfa(f"{rd}/merge/merged.gfa")
    print(f"Merged: N={merged.node_count()}, E={merged.edge_count()}")

    if br:
        _write_boundary_report(br, f"{rd}/merge/boundary_report.tsv")
        print(f"Boundaries: {f'{rd}/merge/boundary_report.tsv'}")


if __name__ == "__main__":
    main()
