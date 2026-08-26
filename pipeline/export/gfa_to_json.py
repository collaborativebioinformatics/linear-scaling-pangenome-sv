"""Export GFA to JSON for the web visualizer."""
import argparse
import json
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph, infer_data_mode


def export(g, label=""):
    return dict(data_mode=infer_data_mode(g), label=label,
                nodes=[dict(id=s.name, len=s.length) for s in g.segments.values()],
                edges=[dict(frm=l.from_node, to=l.to_node) for l in g.links],
                paths={pn: p.segment_names for pn, p in g.paths.items()},
                metrics=dict(nodes=g.node_count(), edges=g.edge_count(),
                             paths=g.path_count(), total_bp=g.total_sequence_bp()))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--output", "-o")
    p.add_argument("--label", "-l", default="")
    a = p.parse_args()
    g = GfaGraph.parse_file(a.input)
    data = export(g, a.label)
    output = a.output or a.input.rsplit(".", 1)[0] + ".json"
    json.dump(data, open(output, "w"), indent=2)
    print(f"{g.node_count()}n -> {output}")


if __name__ == "__main__":
    main()