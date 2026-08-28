"""Graph merging for parallel pangenome graphs.
Strategies: diagnostic_disjoint_union, overlap_aware_stitch."""
from __future__ import annotations
import csv, os, sys
from typing import Dict, List, Optional, Tuple
from pipeline.merge.gfa import (
    GfaGraph, Header, Segment, Link, Path, Walk, _split_orient
)
from pipeline.merge.paths import group_paths_by_haplotype, stitch_haplotype

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
            path=np, tags=dict(w.tags)
        ))


def walks_as_paths(g):
    """Materialise W lines as P lines so a chunk graph has one path view.

    PGGB emits P lines; other builders emit W lines. The stitch works on paths,
    so a W-only chunk graph is converted here instead of being special-cased in
    every downstream step. Names follow PanSN with the subrange the W line
    already carries.
    """
    if g.paths or not g.walks:
        return g
    g = g.copy()
    for w in g.walks:
        name = f"{w.sample}#{w.haplotype}#{w.contig}:{w.start}-{w.end}"
        g.paths[name] = Path(name, list(w.path))
    g.walks = []
    return g


def overlap_aware_stitch(chunk_graphs, ref_name="GRCh38", overlap_bp=100000,
                         chunk_rows=None, chunk_mapping=None):
    """Stitch adjacent chunks into one graph across their reference overlap.

    Each haplotype is cut once inside every overlap and the two sides are
    concatenated, so the overlapping stretch survives exactly once. Nodes are
    keyed by (chunk, segment, sliced offsets), which keeps a segment shared by
    several haplotypes inside a chunk shared after the cut - without that the
    merge would degenerate into one independent chain per haplotype.

    Returns (merged_graph, boundary_reports). A boundary where two chunks do not
    actually overlap on some haplotype is reported FAIL and left unjoined: the
    merge never invents an edge it cannot justify.
    """
    merged = GfaGraph()
    if not chunk_graphs:
        return merged, []

    if chunk_rows is None:
        chunk_rows = _infer_chunk_rows(chunk_graphs)

    chunk_graphs = [(cid, walks_as_paths(g)) for cid, g in chunk_graphs]
    graphs = dict(chunk_graphs)
    merged.headers = [Header(h.version, dict(h.metadata))
                      for h in chunk_graphs[0][1].headers] or [Header("1.1")]

    groups = group_paths_by_haplotype(chunk_graphs, chunk_rows, chunk_mapping)
    stitched, gaps = {}, {}
    for key, entries in groups.items():
        pieces, hap_gaps = stitch_haplotype(chunk_graphs, entries)
        stitched[key] = pieces
        for lc, rc, liv, riv in hap_gaps:
            gaps.setdefault((lc, rc), []).append((key, liv, riv))

    # deterministic node IDs from sorted piece keys, so the output does not
    # depend on the order the chunks were handed to us
    unoriented = sorted({(c, n, lo, hi)
                         for pieces in stitched.values()
                         for c, n, _o, lo, hi in pieces})
    ids = {}
    for i, (cid, name, lo, hi) in enumerate(unoriented, 1):
        nid = f"s{i}"
        ids[(cid, name, lo, hi)] = nid
        merged.segments[nid] = Segment(
            nid, graphs[cid].segments[name].sequence[lo:hi])

    gap_pairs = set(gaps)
    edges = set()

    def _run_steps(run_pieces):
        steps = [ids[(c, n, lo, hi)] + o for c, n, o, lo, hi in run_pieces]
        for a, b in zip(steps, steps[1:]):
            an, ao = _split_orient(a)
            bn, bo = _split_orient(b)
            edges.add((an, ao, bn, bo))
        return steps

    for key in sorted(stitched):
        pieces = stitched[key]
        if not pieces:
            continue
        runs, run = [], [pieces[0]]
        for prev, cur in zip(pieces, pieces[1:]):
            pc, cc = prev[0], cur[0]
            if pc != cc and ((pc, cc) in gap_pairs or (cc, pc) in gap_pairs):
                runs.append(run)
                run = [cur]
            else:
                run.append(cur)
        runs.append(run)
        merged.paths["#".join(key)] = Path("#".join(key), _run_steps(runs[0]))
        for extra in runs[1:]:
            _run_steps(extra)

    # keep chunk-internal links whose endpoints both survived whole: a variant
    # edge that no path in this chunk walks would otherwise be dropped
    whole = {}
    for (cid, name, lo, hi), nid in ids.items():
        if lo == 0 and hi == graphs[cid].segments[name].length:
            whole[(cid, name)] = nid
    for cid, g in chunk_graphs:
        for l in g.links:
            u, v = whole.get((cid, l.from_node)), whole.get((cid, l.to_node))
            if u and v:
                edges.add((u, l.from_orient, v, l.to_orient))

    merged.links = [Link(u, uo, v, vo, "0M") for u, uo, v, vo in sorted(edges)]
    return merged, _boundary_reports(chunk_graphs, groups, gaps, chunk_rows,
                                     overlap_bp)


def _reference_chunk_order(chunk_graphs, groups, rows):
    """Chunk ids sorted by haplotype/reference start, not input order."""
    starts = {}
    for cid, _g in chunk_graphs:
        row = rows.get(cid)
        if row and row.get("reference_start") is not None:
            starts[cid] = int(row["reference_start"])
            continue
        ivs = [iv[0] for entries in groups.values()
               for c, _p, iv in entries if c == cid]
        starts[cid] = min(ivs) if ivs else 0
    return sorted((cid for cid, _ in chunk_graphs),
                  key=lambda c: (starts[c], c))


def _failed_haplotypes(left, right, order, gaps):
    """Haplotypes that fail this boundary, including a gap that spans it.

    stitch_haplotype keys gaps by the two surviving windows, which may skip a
    dropped middle chunk. A consecutive pair inside that span must still FAIL.
    """
    try:
        i, j = order.index(left), order.index(right)
    except ValueError:
        return list(gaps.get((left, right), []))
    if i > j:
        i, j = j, i
    seen = {}
    for (gl, gr), ghaps in gaps.items():
        try:
            gi, gj = order.index(gl), order.index(gr)
        except ValueError:
            continue
        if gi > gj:
            gi, gj = gj, gi
        if gi <= i and j <= gj:
            for item in ghaps:
                seen[item[0]] = item
    return list(seen.values())


def _boundary_reports(chunk_graphs, groups, gaps, chunk_rows, overlap_bp):
    rows = chunk_rows or {}
    order = _reference_chunk_order(chunk_graphs, groups, rows)
    reports = []
    for left, right in zip(order, order[1:]):
        actual = _actual_overlap(rows.get(left), rows.get(right), overlap_bp)
        shared = sorted(k for k, e in groups.items()
                        if {left, right} <= {cid for cid, _p, _i in e})
        failed = _failed_haplotypes(left, right, order, gaps)
        failed_keys = {g[0] for g in failed}
        joined = [k for k in shared if k not in failed_keys]
        reports.append({
            "boundary": f"{left}--{right}",
            "left_chunk": left,
            "right_chunk": right,
            "reference_overlap_bp": actual,
            "anchor_found": bool(joined),
            "haplotypes_preserved": not failed,
            "haplotypes_joined": len(joined),
            "haplotypes_unjoined": len(failed),
            "status": "PASS" if joined and not failed else
                      ("FAIL" if failed else "WARN"),
            "message": _boundary_message(joined, failed),
        })
    return reports


def _boundary_message(joined, failed):
    if failed:
        names = ", ".join("#".join(k) for k, _l, _r in failed[:3])
        return f"no overlap for {len(failed)} haplotype(s): {names}"
    if not joined:
        return "no haplotype present in both chunks"
    return f"stitched {len(joined)} haplotype(s) at overlap midpoint"


def _actual_overlap(left_row, right_row, fallback):
    """Overlap two chunk windows really share on the reference.

    make_chunks treats overlap_bp as the pairwise overlap: each interior
    core is padded by overlap/2 per side. The manifest windows are the
    source of truth, so a config fallback is only used when rows are missing.
    """
    if not left_row or not right_row:
        return fallback
    lo = max(int(left_row["reference_start"]), int(right_row["reference_start"]))
    hi = min(int(left_row["reference_end"]), int(right_row["reference_end"]))
    return max(0, hi - lo)


def _load_chunks(cm_path, chunk_dir=None):
    """([(chunk_id, GfaGraph)], {chunk_id: manifest_row}) for chunks already built.

    chunk_dir defaults to a `chunks/` folder next to the manifest, otherwise
    the manifest's own directory (covers both work/chunks/ and work/demo/).
    """
    if not os.path.exists(cm_path):
        return [], {}
    if chunk_dir is None:
        d = os.path.dirname(os.path.abspath(cm_path)) or "."
        nested = os.path.join(d, "chunks")
        chunk_dir = nested if os.path.isdir(nested) else d
    result, rows = [], {}
    with open(cm_path) as f:
        for row in csv.DictReader(f, delimiter=T):
            gp = os.path.join(chunk_dir, f"{row['chunk_id']}.gfa")
            if os.path.exists(gp):
                result.append((row["chunk_id"], GfaGraph.parse_file(gp)))
                rows[row["chunk_id"]] = row
    return result, rows


def _infer_chunk_rows(chunk_graphs):
    """Load a nearby chunk_manifest.tsv when the caller omitted chunk_rows."""
    for _cid, g in chunk_graphs:
        src = getattr(g, "source", None)
        if not src:
            continue
        d = os.path.dirname(os.path.abspath(src))
        for cand in (os.path.join(d, "chunk_manifest.tsv"),
                     os.path.join(os.path.dirname(d), "chunk_manifest.tsv")):
            if os.path.exists(cand):
                _loaded, rows = _load_chunks(cand, chunk_dir=d)
                if rows:
                    return rows
    return {}


def _load_chunk_mapping(path="results/preparation/chunk_mapping.tsv"):
    """Rows from chunk_mapping.tsv (source_start/source_end/strand per hap)."""
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter=T))


def _write_boundary_report(boundaries, path):
    keys = [
        "boundary", "left_chunk", "right_chunk", "reference_overlap_bp",
        "anchor_found", "haplotypes_preserved", "haplotypes_joined",
        "haplotypes_unjoined", "status", "message"
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

    chunks, rows = _load_chunks("work/chunks/chunk_manifest.tsv")
    if not chunks:
        print("No chunks. Run make chunks.")
        return
    print(f"Loaded {len(chunks)} chunks")

    if strategy == "disjoint_union":
        merged = diagnostic_disjoint_union(chunks)
        br = []
    else:
        obp = config.get("parallel", {}).get("overlap_bp", 100000)
        mapping = _load_chunk_mapping()
        merged, br = overlap_aware_stitch(
            chunks, overlap_bp=obp, chunk_rows=rows, chunk_mapping=mapping)

    os.makedirs(f"{rd}/merge", exist_ok=True)
    merged.write_gfa(f"{rd}/merge/merged.gfa")
    print(f"Merged: N={merged.node_count()}, E={merged.edge_count()}")

    if br:
        _write_boundary_report(br, f"{rd}/merge/boundary_report.tsv")
        failed = [b for b in br if b["status"] == "FAIL"]
        print(f"Boundaries: {rd}/merge/boundary_report.tsv "
              f"({len(br) - len(failed)}/{len(br)} PASS)")
        if failed and config.get("merge", {}).get("fail_on_unresolved_boundary"):
            sys.exit(1)


if __name__ == "__main__":
    main()
