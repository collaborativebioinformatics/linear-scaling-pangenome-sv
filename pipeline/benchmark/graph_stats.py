"""Compute graph statistics for baseline and merged GFAs."""
import json
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph


def stats(g, label=""):
    return dict(label=label, nodes=g.node_count(), edges=g.edge_count(),
                paths=g.path_count(), walks=g.walk_count(),
                total_bp=g.total_sequence_bp(),
                samples=sorted(g.get_sample_names()))


def main():
    rd = "results"
    targets = [(f"{rd}/baseline/baseline.gfa", "baseline"),
               (f"{rd}/merge/merged.gfa", "merged")]
    out = {}
    for path, label in targets:
        if os.path.exists(path):
            g = GfaGraph.parse_file(path)
            out[label] = stats(g, label)
            print(f"{label}: {g.node_count()}n {g.edge_count()}e")
    if out:
        op = f"{rd}/benchmark/graph_metrics.json"
        os.makedirs(os.path.dirname(op), exist_ok=True)
        json.dump(out, open(op, "w"), indent=2)
        print(f"Metrics: {op}")


if __name__ == "__main__":
    main()