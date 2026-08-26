"""
paths.py — Path traversal across chunk boundaries.

A chunk graph covers one reference window; adjacent windows overlap on
purpose. Every haplotype therefore appears once per chunk it touches, and the
overlapping stretch is present twice. Stitching one haplotype back together is
three steps:

  1. locate each chunk-local path on the haplotype's own coordinates
     (`path_interval`)
  2. pick a single cut point inside each overlap (`seam_position`)
  3. keep only the part of each chunk path on its side of the cut
     (`slice_steps`)

Slicing is done in *segment offsets*, not by minting per-haplotype sequence, so
two haplotypes that share a node inside a chunk still share it after the cut.
That is what keeps the merged graph a pangenome instead of a bundle of
independent chains.
"""
from pipeline.merge.gfa import GfaGraph, haplotype_key, parse_pansn


def path_interval(graph, path_name, chunk_row=None):
    """(start, end) of a chunk-local path on its own contig.

    Prefers the `:start-end` subrange PGGB carries over from a sliced input
    FASTA. Falls back to the chunk manifest's reference window, which is only
    correct when haplotype coordinates equal reference coordinates (synthetic
    demo data). Length always comes from the spelled path, never from the
    manifest, so a path that stops short reports where it really stops.
    """
    _s, _h, _c, start, end, _cid = parse_pansn(path_name)
    if start is None:
        if chunk_row is None:
            raise ValueError(
                f"path {path_name} has no subrange and no chunk manifest row")
        start = int(chunk_row["reference_start"])
    return start, start + graph.path_length(path_name)


def seam_position(left, right):
    """Cut point inside the overlap of two adjacent intervals, or None.

    None means the two chunks do not actually overlap on this haplotype: a
    boundary the merge cannot resolve, which the caller must report rather
    than paper over.
    """
    lo, hi = max(left[0], right[0]), min(left[1], right[1])
    if hi <= lo:
        return None
    return (lo + hi) // 2


def slice_steps(graph, steps, start, keep_from, keep_to):
    """Trim an oriented step list to the coordinate window [keep_from, keep_to).

    `steps` is [(segment_name, orient)], `start` is the coordinate of its first
    base. Returns [(segment_name, orient, off_lo, off_hi)] where the offsets are
    into the segment's *forward* sequence, so a reverse step keeps the piece
    that actually falls in the window.
    """
    out, pos = [], start
    for name, orient in steps:
        length = graph.segments[name].length
        lo, hi = max(pos, keep_from), min(pos + length, keep_to)
        if hi > lo:
            rel_lo, rel_hi = lo - pos, hi - pos
            if orient == "+":
                out.append((name, orient, rel_lo, rel_hi))
            else:
                out.append((name, orient, length - rel_hi, length - rel_lo))
        pos += length
    return out


def group_paths_by_haplotype(chunk_graphs, chunk_rows=None):
    """{(sample, hap, contig): [(chunk_id, path_name, (start, end))]} in order.

    chunk_graphs is [(chunk_id, GfaGraph)]; chunk_rows maps chunk_id to its
    manifest row and may be omitted when the path names carry subranges.
    """
    rows = chunk_rows or {}
    groups = {}
    for cid, g in chunk_graphs:
        for pn in g.paths:
            iv = path_interval(g, pn, rows.get(cid))
            groups.setdefault(haplotype_key(pn), []).append((cid, pn, iv))
    for key in groups:
        groups[key].sort(key=lambda t: (t[2][0], t[0]))
    return groups


def stitch_haplotype(chunk_graphs, entries):
    """Cut one haplotype's chunk paths at their seams and concatenate.

    Returns (pieces, gaps) where pieces is
    [(chunk_id, segment_name, orient, off_lo, off_hi)] in haplotype order and
    gaps is [(left_chunk, right_chunk, left_interval, right_interval)] for each
    pair that failed to overlap, including surviving windows left non-adjacent
    by a skipped (empty keep) chunk.
    """
    graphs = dict(chunk_graphs)
    seams, gaps = [], []
    for (lc, _lp, liv), (rc, _rp, riv) in zip(entries, entries[1:]):
        seam = seam_position(liv, riv)
        if seam is None:
            gaps.append((lc, rc, liv, riv))
        seams.append(seam)

    pieces = []
    last_keep_to = last_cid = last_iv = None
    known_gaps = {(g[0], g[1]) for g in gaps}
    for i, (cid, pn, (start, end)) in enumerate(entries):
        left_seam = seams[i - 1] if i > 0 else None
        right_seam = seams[i] if i < len(seams) else None
        keep_from = start if left_seam is None else max(start, left_seam)
        keep_to = end if right_seam is None else min(end, right_seam)
        if keep_to <= keep_from:
            continue
        # A skipped (empty keep) chunk can leave two surviving windows that
        # do not meet. Concatenating those pieces would invent an edge.
        if last_keep_to is not None and keep_from != last_keep_to:
            pair = (last_cid, cid)
            if pair not in known_gaps:
                gaps.append((last_cid, cid, last_iv, (start, end)))
                known_gaps.add(pair)
        g = graphs[cid]
        for name, orient, lo, hi in slice_steps(
                g, g.path_steps(pn), start, keep_from, keep_to):
            pieces.append((cid, name, orient, lo, hi))
        last_keep_to, last_cid, last_iv = keep_to, cid, (start, end)
    return pieces, gaps


def compare_path_lengths(baseline: GfaGraph, merged: GfaGraph,
                         path_name: str) -> dict:
    bp = baseline.paths.get(path_name)
    mp = merged.paths.get(path_name)
    return {
        "path": path_name,
        "baseline_segments": len(bp.segment_names) if bp else 0,
        "merged_segments": len(mp.segment_names) if mp else 0,
        "in_baseline": bp is not None,
        "in_merged": mp is not None,
    }
