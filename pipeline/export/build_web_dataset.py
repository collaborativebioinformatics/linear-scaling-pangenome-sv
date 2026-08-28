"""build_web_dataset.py — Convert GFA results into bounded compact JSON for the web.

This is the ONLY contract between genomics results and the browser frontend.

Inputs (all optional — whatever exists is discovered dynamically):
    results/baseline/baseline.gfa     (monolithic baseline)
    results/merge/merged.gfa          (stitched merged)
    work/demo/baseline.gfa            (synthetic fallback)
    work/demo/merged.gfa              (synthetic fallback)
    work/chunks/chunk_manifest.tsv    (chunk coordinate windows)
    work/demo/chunk_manifest.tsv      (synthetic fallback)

Outputs:
    results/web/manifest.json                       sample/graph discovery + status
    results/web/overview.json                       chunk coordinate windows
    results/web/graphs/<graph>/<sample>_<hap>.json  bounded per-haplotype graph

Safety limits (env-overridable):
    WEB_GRAPH_MAX_NODES       default 1500
    WEB_GRAPH_MAX_EDGES       default 4000
    WEB_GRAPH_MAX_PATH_STEPS  default 5000

No raw nucleotide sequence is ever written to a web JSON file.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph, parse_pansn, _split_orient

WEB_GRAPH_MAX_NODES = int(os.environ.get("WEB_GRAPH_MAX_NODES", "1500"))
WEB_GRAPH_MAX_EDGES = int(os.environ.get("WEB_GRAPH_MAX_EDGES", "4000"))
WEB_GRAPH_MAX_PATH_STEPS = int(os.environ.get("WEB_GRAPH_MAX_PATH_STEPS", "5000"))


def _samples_from_graph(g):
    """Discover {sample: set(haplotypes)} from P and W records."""
    out = {}
    for pn in g.paths:
        sample, hap, contig, _s, _e, _c = parse_pansn(pn)
        out.setdefault(sample, set()).add(hap)
    for w in g.walks:
        out.setdefault(w.sample, set()).add(w.haplotype)
    return out


def _build_selected_path(g, sample, hap):
    """Return (path_name, [(seg, orient)]) for the best-matching path/walk."""
    for pn in g.paths:
        s, h, contig, _a, _b, _c = parse_pansn(pn)
        if s == sample and h == hap:
            return pn, g.path_steps(pn)
    for w in g.walks:
        if w.sample == sample and w.haplotype == hap:
            return f"{w.sample}#{w.haplotype}#{w.contig}", g.walk_steps(w)
    return None, []


def _extract_sample_graph(g, sample, hap, graph_label):
    """Build a bounded per-haplotype graph: selected path + one-hop neighbors."""
    path_name, steps = _build_selected_path(g, sample, hap)
    selected = set()
    for seg, orient in steps:
        selected.add(seg)

    include = set(selected)
    for link in g.links:
        if link.from_node in selected:
            include.add(link.to_node)
        if link.to_node in selected:
            include.add(link.from_node)

    degree = {}
    neighbors = {}
    for link in g.links:
        if link.from_node in include or link.to_node in include:
            degree[link.from_node] = degree.get(link.from_node, 0) + 1
            degree[link.to_node] = degree.get(link.to_node, 0) + 1
            neighbors.setdefault(link.from_node, set()).add(link.to_node)
            neighbors.setdefault(link.to_node, set()).add(link.from_node)

    is_reference = (sample == "GRCh38")
    nodes = []
    for seg_name in include:
        seg = g.segments.get(seg_name)
        if seg is None:
            continue
        nodes.append({
            "id": seg_name,
            "length": seg.length,
            "on_selected_path": seg_name in selected,
            "on_reference": is_reference and (seg_name in selected),
            "degree": degree.get(seg_name, 0),
            "neighbors": sorted(neighbors.get(seg_name, [])),
        })

    edges = []
    seen_edges = set()
    for link in g.links:
        if link.from_node not in include or link.to_node not in include:
            continue
        key = (link.from_node, link.from_orient, link.to_node, link.to_orient)
        if key in seen_edges:
            continue
        seen_edges.add(key)
        edges.append({
            "source": link.from_node,
            "target": link.to_node,
            "source_orientation": link.from_orient,
            "target_orientation": link.to_orient,
            "on_selected_path": link.from_node in selected and link.to_node in selected,
            "on_reference": is_reference and link.from_node in selected and link.to_node in selected,
        })

    original_nodes = len(nodes)
    original_edges = len(edges)
    truncated = original_nodes > WEB_GRAPH_MAX_NODES or original_edges > WEB_GRAPH_MAX_EDGES
    nodes = nodes[:WEB_GRAPH_MAX_NODES]
    edges = edges[:WEB_GRAPH_MAX_EDGES]

    path_steps = [{"node": seg, "orientation": orient} for seg, orient in steps]
    path_steps = path_steps[:WEB_GRAPH_MAX_PATH_STEPS]
    path_len_bp = sum(g.segments[s].length for s, _ in steps if s in g.segments)

    return {
        "schema_version": "1",
        "graph": graph_label,
        "sample": sample,
        "haplotype": hap,
        "path_name": path_name,
        "nodes": nodes,
        "edges": edges,
        "path": {"steps": path_steps, "length_bp": path_len_bp},
        "metrics": {"nodes": len(nodes), "edges": len(edges),
                    "path_steps": len(path_steps)},
        "truncated": truncated,
        "original_counts": {"nodes": original_nodes, "edges": original_edges,
                            "path_steps": len(steps)},
    }


def _discover_graphs():
    """Locate baseline + merged GFAs across real and synthetic locations."""
    candidates = {
        "baseline": ["results/baseline/baseline.gfa", "work/demo/baseline.gfa"],
        "merged": ["results/merge/merged.gfa", "work/demo/merged.gfa"],
    }
    found = {}
    for label, paths in candidates.items():
        for p in paths:
            if os.path.exists(p):
                found[label] = p
                break
    return found


def _discover_chunk_manifest():
    for p in ["work/chunks/chunk_manifest.tsv", "work/demo/chunk_manifest.tsv"]:
        if os.path.exists(p):
            return p
    return None


def _load_chunk_rows(manifest_path):
    rows = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    return rows


def _manifest_samples(samples_by_graph):
    merged = {}
    for graph_label, samples in samples_by_graph.items():
        for sample, haps in samples.items():
            merged.setdefault(sample, set()).update(haps)
    return [{"sample": s, "haplotypes": sorted(hs)}
            for s, hs in sorted(merged.items())]


def _hap_label(hap):
    return {"0": "reference", "1": "paternal", "2": "maternal"}.get(hap, hap)


def main():
    os.makedirs("results/web", exist_ok=True)
    graphs = _discover_graphs()
    chunk_manifest = _discover_chunk_manifest()

    samples_by_graph = {}
    graph_meta = {}
    for label, path in graphs.items():
        g = GfaGraph.parse_file(path)
        samples_by_graph[label] = _samples_from_graph(g)
        graph_meta[label] = {
            "nodes": g.node_count(),
            "edges": g.edge_count(),
            "paths": g.path_count(),
            "walks": g.walk_count(),
        }

    samples = _manifest_samples(samples_by_graph)
    has_real_baseline = os.path.exists("results/baseline/baseline.gfa")
    data_mode = "real" if has_real_baseline else "synthetic"

    os.makedirs("results/web/graphs", exist_ok=True)
    written = 0
    for label, path in graphs.items():
        g = GfaGraph.parse_file(path)
        d = os.path.join("results/web/graphs", label)
        os.makedirs(d, exist_ok=True)
        for sample, haps in samples_by_graph[label].items():
            for hap in haps:
                sg = _extract_sample_graph(g, sample, hap, label)
                fn = os.path.join(d, f"{sample}_{hap}.json")
                with open(fn, "w") as f:
                    json.dump(sg, f)
                written += 1

    chunks = []
    if chunk_manifest:
        for r in _load_chunk_rows(chunk_manifest):
            chunks.append({
                "chunk_id": r.get("chunk_id", ""),
                "chrom": r.get("chrom", ""),
                "reference_start": int(r.get("reference_start", 0) or 0),
                "reference_end": int(r.get("reference_end", 0) or 0),
                "core_start": int(r.get("core_start", 0) or 0),
                "core_end": int(r.get("core_end", 0) or 0),
                "overlap_left": int(r.get("overlap_left", 0) or 0),
                "overlap_right": int(r.get("overlap_right", 0) or 0),
                "pairwise_overlap_bp": int(r.get("pairwise_overlap_bp", 0) or 0),
            })

    pipeline_status = {
        "baseline": "IMPLEMENTED" if "baseline" in graphs else "NOT_AVAILABLE",
        "parallel_chunks": "IMPLEMENTED" if chunks else "NOT_AVAILABLE",
        "stitch": "IMPLEMENTED" if "merged" in graphs else "NOT_IMPLEMENTED",
        "equivalence": "NOT_RUN",
    }

    manifest = {
        "schema_version": "1",
        "run_id": "latest",
        "data_mode": data_mode,
        "target": {"reference": "GRCh38", "chromosome": "chr21",
                   "start": 0, "end": 0},
        "pipeline_status": pipeline_status,
        "samples": [{"sample": s["sample"],
                     "haplotypes": s["haplotypes"],
                     "hap_labels": {h: _hap_label(h) for h in s["haplotypes"]}}
                    for s in samples],
        "graphs": graph_meta,
    }
    with open("results/web/manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    overview = {
        "schema_version": "1",
        "target": manifest["target"],
        "chunks": chunks,
        "graphs": graph_meta,
        "pipeline_status": pipeline_status,
    }
    with open("results/web/overview.json", "w") as f:
        json.dump(overview, f, indent=2)

    print(f"build_web_dataset: {len(graphs)} graph(s), {len(samples)} sample(s), "
          f"{written} sample graph JSON(s), {len(chunks)} chunk(s)")
    print("  manifest -> results/web/manifest.json")
    print("  overview -> results/web/overview.json")


if __name__ == "__main__":
    main()
