"""select_visualization_region.py — Pick a fixed, informative subregion for Cytoscape.

Scans a GFA for a bounded, branch-dense window suitable for visualization.
The FULL benchmark metrics stay on the whole region; this only chooses the
small window the browser renders.

Deterministic: fixed window size, fixed step, lexicographic tie-break.

Usage:
    python3 pipeline/export/select_visualization_region.py <gfa> \
        [--window-bp 50000] [--step-bp 10000] [--out results/web/visualization_region.json]
"""
import argparse
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph, parse_pansn


def _node_sort_key(name):
    try:
        return (0, int(name))
    except ValueError:
        return (1, name)


def _path_node_sets(g):
    """{sample#hap: set(node_ids)} for every path/walk."""
    out = {}
    for pn in g.paths:
        s, h, c, _a, _b, _c2 = parse_pansn(pn)
        out.setdefault(f"{s}#{h}", set()).update(
            n for n, _ in g.path_steps(pn)
        )
    for w in g.walks:
        key = f"{w.sample}#{w.haplotype}"
        out.setdefault(key, set()).update(n for n, _ in g.walk_steps(w))
    return out


def score_window(g, node_ids, path_sets, ref_nodes):
    """Score a set of nodes for visual interest.

    branch density + path disagreement + sample diversity, penalizing
    unrenderable node counts.
    """
    n = len(node_ids)
    if n == 0:
        return 0.0

    degree = defaultdict(int)
    for link in g.links:
        if link.from_node in node_ids:
            degree[link.from_node] += 1
        if link.to_node in node_ids:
            degree[link.to_node] += 1

    branch_nodes = sum(1 for d in degree.values() if d > 2)
    branch_density = branch_nodes / n

    # sample diversity: how many distinct haplotypes traverse this window
    samples_present = set()
    for key, nodes in path_sets.items():
        if nodes & node_ids:
            samples_present.add(key)
    sample_diversity = len(samples_present) / max(1, len(path_sets))

    # path disagreement: nodes used by exactly one haplotype (unique)
    usage = defaultdict(int)
    for key, nodes in path_sets.items():
        for nd in (nodes & node_ids):
            usage[nd] += 1
    unique = sum(1 for c in usage.values() if c == 1)
    path_diversity = unique / n

    ref_fraction = len(node_ids & ref_nodes) / n

    # penalize huge node counts (unrenderable) and low total bp
    total_bp = sum(g.segments[nd].length for nd in node_ids if nd in g.segments)
    size_penalty = min(1.0, n / 800) * 0.3

    score = (
        3.0 * branch_density
        + 3.0 * path_diversity
        + 2.0 * sample_diversity
        + 1.0 * min(1.0, ref_fraction + 0.3)
        - size_penalty
    )
    return score


def select_region(gfa_path, window_bp=50000, step_bp=10000):
    g = GfaGraph.parse_file(gfa_path)
    path_sets = _path_node_sets(g)
    ref_nodes = set()
    for pn in g.paths:
        s, h, c, _a, _b, _c2 = parse_pansn(pn)
        if s == "GRCh38":
            ref_nodes.update(n for n, _ in g.path_steps(pn))

    all_nodes = sorted(g.segments.keys(), key=_node_sort_key)
    n = len(all_nodes)

    # Sliding window over node index space (nodes are roughly in genomic order
    # for a linear reference backbone, but we score on graph structure so it
    # works for any ordering).
    windows = []
    start = 0
    while start < n:
        win = set(all_nodes[start:start + max(1, window_bp // 500)])
        # Better: use bp-based window over cumulative path bp when possible,
        # but node-based windows are deterministic and graph-aware.
        if win:
            sc = score_window(g, win, path_sets, ref_nodes)
            branch_nodes = sum(
                1 for link in g.links
                if link.from_node in win and link.to_node in win
                and (sum(1 for l in g.links if l.from_node == link.from_node) > 2)
            )
            samples_present = sum(
                1 for nodes in path_sets.values() if nodes & win
            )
            windows.append({
                "node_start": start,
                "node_end": min(start + max(1, window_bp // 500), n),
                "node_count": len(win),
                "score": round(sc, 4),
                "branch_nodes": branch_nodes,
                "samples_present": samples_present,
            })
        start += max(1, step_bp // 500)

    if not windows:
        return None

    windows.sort(key=lambda w: (-w["score"], w["node_start"]))
    best = windows[0]

    # Compute bp span using segments in the best window
    win_nodes = set(all_nodes[best["node_start"]:best["node_end"]])
    total_bp = sum(g.segments[nd].length for nd in win_nodes if nd in g.segments)

    return {
        "gfa_source": gfa_path,
        "window_bp": window_bp,
        "step_bp": step_bp,
        "node_start": best["node_start"],
        "node_end": best["node_end"],
        "node_count": best["node_count"],
        "score": best["score"],
        "branch_nodes": best["branch_nodes"],
        "samples_present": best["samples_present"],
        "approx_bp": total_bp,
        "selection_method":
            "highest branch/path-diversity score within benchmark graph",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gfa")
    ap.add_argument("--window-bp", type=int, default=50000)
    ap.add_argument("--step-bp", type=int, default=10000)
    ap.add_argument("--out", default="results/web/visualization_region.json")
    args = ap.parse_args()

    result = select_region(args.gfa, args.window_bp, args.step_bp)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
