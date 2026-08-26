"""
find_bubbles.py — First-pass bubble detection in GFA.
TODO (Quang/Michael): Complex nested bubble detection.
TODO (Ali): Integrate with web display.
"""
import argparse
import json
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph


def find_bubbles(graph: GfaGraph, max_depth: int = 5) -> list:
    out = {}
    for link in graph.links:
        out.setdefault(link.from_node, []).append(link.to_node)

    def reachable(start: str, depth: int) -> set:
        if depth <= 0:
            return {start}
        r = {start}
        for n in out.get(start, []):
            r |= reachable(n, depth - 1)
        return r

    bubbles = []
    for node, targets in out.items():
        if len(targets) >= 2:
            for i, t1 in enumerate(targets):
                for t2 in targets[i + 1:]:
                    common = reachable(t1, max_depth) & reachable(t2, max_depth)
                    if common:
                        sink = min(common)
                        bubbles.append({
                            "source": node,
                            "sink": sink,
                            "alleles": [t1, t2],
                            "depth": max_depth,
                        })
    return bubbles


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--output", "-o")
    args = parser.parse_args()
    graph = GfaGraph.parse_file(args.input)
    bubbles = find_bubbles(graph)
    output = args.output or args.input.rsplit(".", 1)[0] + "_bubbles.json"
    with open(output, "w") as f:
        json.dump(bubbles, f, indent=2)
    print(f"Found {len(bubbles)} bubbles -> {output}")


if __name__ == "__main__":
    main()