"""Phase 1-2: Diagnose zero-edge chunk GFAs and compare with baseline."""
import sys, json
from collections import Counter
sys.path.insert(0, ".")
from pipeline.merge.gfa import GfaGraph

print("=" * 60)
print("PHASE 1 — CHUNK GFA DIAGNOSTICS")
print("=" * 60)

diag = []

for i in [1, 2, 3]:
    path = f"results/chunks/chunk_000{i}.gfa"
    g = GfaGraph.parse_file(path)
    lines = open(path).readlines()
    counts = Counter(l[0] for l in lines if l.strip())

    info = {
        "chunk_id": f"chunk_000{i}",
        "file_lines": len(lines),
        "record_counts": dict(counts),
        "segments": g.node_count(),
        "links": g.edge_count(),
        "paths": g.path_count(),
        "walks": g.walk_count(),
        "path_names": sorted(g.paths.keys()),
        "total_segment_bp": sum(s.length for s in g.segments.values()),
        "status": "INVALID_ZERO_EDGE_GRAPH" if g.edge_count() == 0 else "VALID",
    }
    diag.append(info)

    print(f"\n=== chunk_000{i} ===")
    print(f"  Lines: {len(lines)}, Records: {dict(counts)}")
    print(f"  Segments: {g.node_count()}, Links(L/E): {g.edge_count()}, Paths: {g.path_count()}")
    print(f"  Total bp: {info['total_segment_bp']:,}")
    print(f"  Status: {info['status']}")
    for pn in sorted(g.paths.keys()):
        steps = g.path_steps(pn)
        bp = sum(g.segments[n].length for n, _ in steps if n in g.segments)
        print(f"    {pn}: {len(steps)} steps, {bp} bp")

print("\n" + "=" * 60)
print("BASELINE COMPARISON")
print("=" * 60)

bg = GfaGraph.parse_file("results/baseline/baseline.gfa")
blines = open("results/baseline/baseline.gfa").readlines()
bcounts = Counter(l[0] for l in blines if l.strip())
print(f"  Lines: {len(blines)}, Records: {dict(bcounts)}")
print(f"  Segments: {bg.node_count()}, Links(L): {bg.edge_count()}, Paths: {bg.path_count()}")
for pn in sorted(bg.paths.keys()):
    steps = bg.path_steps(pn)
    bp = sum(bg.segments[n].length for n, _ in steps if n in bg.segments)
    print(f"    {pn}: {len(steps)} steps, {bp:,} bp")

# Save diagnostics
import os
os.makedirs("results/validation", exist_ok=True)
with open("results/validation/chunk_graph_diagnostics.json", "w") as f:
    json.dump(diag, f, indent=2)
print("\nSaved: results/validation/chunk_graph_diagnostics.json")