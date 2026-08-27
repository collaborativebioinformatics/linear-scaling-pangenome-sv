"""Export GFA to JSON for the web visualizer.
Bounded export: large graphs are truncated to prevent memory issues.
"""
import argparse
import json
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph, infer_data_mode


def export(g, label="", max_nodes=5000, max_edges=15000, max_walks=50):
    """Export graph to JSON with bounded size limits."""
    nodes = [dict(id=s.name, len=s.length) for s in g.segments.values()]
    edges = [dict(frm=l.from_node, to=l.to_node) for l in g.links]
    paths = {pn: p.segment_names for pn, p in g.paths.items()}
    walks = [dict(sample=w.sample, haplotype=w.haplotype,
                  contig=w.contig, start=w.start, end=w.end,
                  path=w.path) for w in g.walks]

    # Truncate if over limits
    truncated = False
    if len(nodes) > max_nodes:
        nodes = nodes[:max_nodes]
        truncated = True
    if len(edges) > max_edges:
        edges = edges[:max_edges]
        truncated = True
    if len(walks) > max_walks:
        walks = walks[:max_walks]
        truncated = True

    return dict(data_mode=infer_data_mode(g), label=label,
                nodes=nodes, edges=edges,
                paths=paths, walks=walks,
                truncated=truncated,
                metrics=dict(nodes=g.node_count(), edges=g.edge_count(),
                             paths=g.path_count(), walks=g.walk_count(),
                             total_bp=g.total_sequence_bp()))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--output", "-o")
    p.add_argument("--label", "-l", default="")
    p.add_argument("--max-nodes", type=int, default=5000)
    p.add_argument("--max-edges", type=int, default=15000)
    a = p.parse_args()
    g = GfaGraph.parse_file(a.input)
    data = export(g, a.label, a.max_nodes, a.max_edges)
    output = a.output or a.input.rsplit(".", 1)[0] + ".json"
    json.dump(data, open(output, "w"), indent=2)
    status = "truncated" if data["truncated"] else "full"
    print(f"{g.node_count()}n, {g.edge_count()}e, {g.walk_count()}w -> {output} ({status})")


if __name__ == "__main__":
    main()