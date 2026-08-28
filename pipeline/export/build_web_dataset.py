"""build_web_dataset.py — Convert GFA results into bounded compact JSON for the web.

This is the ONLY contract between genomics results and the browser frontend.

Schema v2:
  - on_reference -> is_shared (multi-sample intersection)
  - Adds on_reference_path (actual GRCh38 path membership)
  - Adds sample_count (distinct sample-haplotypes traversing node)
  - Graph-level overview.json (stable shared topology)
  - Per-haplotype focus JSON (selected path + one-hop)
  - Path diagnostics + provenance recording

Inputs (all optional — discovered dynamically):
    results/baseline/baseline.gfa     (monolithic baseline)
    results/merge/merged.gfa          (stitched graph)
    work/demo/baseline.gfa            (synthetic fallback)
    work/demo/merged.gfa              (synthetic fallback)
    work/chunks/chunk_manifest.tsv    (chunk coordinate windows)

Safety limits (env-overridable):
    WEB_GRAPH_MAX_NODES       default 4000
    WEB_GRAPH_MAX_EDGES       default 10000
    WEB_GRAPH_MAX_PATH_STEPS  default 20000

No raw nucleotide sequence is ever written to a web JSON file.
"""
import csv
import hashlib
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph, parse_pansn

WEB_GRAPH_MAX_NODES = int(os.environ.get("WEB_GRAPH_MAX_NODES", "4000"))
WEB_GRAPH_MAX_EDGES = int(os.environ.get("WEB_GRAPH_MAX_EDGES", "10000"))
WEB_GRAPH_MAX_PATH_STEPS = int(os.environ.get("WEB_GRAPH_MAX_PATH_STEPS", "20000"))


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


def _build_shared_nodes(g):
    """Return set of node IDs that appear on paths from multiple samples.
    
    In PGGB graphs, the reference may be a single isolated node (node 1).
    Instead of tying "reference" to that node, we identify shared/common
    nodes — nodes that appear on paths from at least two different individuals.
    """
    sample_nodes = {}  # sample -> set of nodes
    for pn in g.paths:
        s, h, c, _a, _b, _c2 = parse_pansn(pn)
        nodes = set()
        common_set = sample_nodes.setdefault(s, set())
        for seg, orient in g.path_steps(pn):
            nodes.add(seg)
        common_set.update(nodes)
    for w in g.walks:
        nodes = set()
        common_set = sample_nodes.setdefault(w.sample, set())
        for seg, orient in g.walk_steps(w):
            nodes.add(seg)
        common_set.update(nodes)
    
    # Nodes that appear in >1 sample
    shared = set()
    samples_list = list(sample_nodes.values())
    for i in range(len(samples_list)):
        for j in range(i + 1, len(samples_list)):
            shared |= samples_list[i] & samples_list[j]
    return shared


def _build_ref_nodes(g):
    """Nodes on the actual GRCh38 reference path."""
    ref_nodes = set()
    for pn in g.paths:
        sample, hap, contig, _s, _e, _c = parse_pansn(pn)
        if sample == "GRCh38":
            for seg, orient in g.path_steps(pn):
                ref_nodes.add(seg)
    for w in g.walks:
        if w.sample == "GRCh38":
            for seg, orient in g.walk_steps(w):
                ref_nodes.add(seg)
    return ref_nodes


def _sample_count_map(g):
    """{node_id: count of distinct sample-haplotypes traversing it}."""
    node_samples = {}
    for pn in g.paths:
        s, h, c, _a, _b, _c2 = parse_pansn(pn)
        key = f"{s}#{h}"
        for seg, orient in g.path_steps(pn):
            node_samples.setdefault(seg, set()).add(key)
    for w in g.walks:
        key = f"{w.sample}#{w.haplotype}"
        for seg, orient in g.walk_steps(w):
            node_samples.setdefault(seg, set()).add(key)
    return {node: len(ss) for node, ss in node_samples.items()}


def _node_sort_key(name):
    """Sort nodes deterministically: try integer, else lexicographic."""
    try:
        return (0, int(name))
    except ValueError:
        return (1, name)


def _git_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()[:12]
    except Exception:
        return "unknown"


def _file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_config_target():
    try:
        import yaml
        with open("config/pipeline.yaml") as f:
            cfg = yaml.safe_load(f)
        t = cfg.get("target", {})
        return {
            "reference": t.get("reference", "GRCh38"),
            "chromosome": t.get("chromosome", "chr21"),
            "start": int(t.get("start", 0)),
            "end": int(t.get("end", 0)),
        }
    except Exception:
        return {"reference": "GRCh38", "chromosome": "chr21", "start": 0, "end": 0}


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

    # Compute shared nodes from the full graph for coloring only.
    shared_path = _build_shared_nodes(g)

    degree = {}
    neighbors = {}
    for link in g.links:
        if link.from_node in include or link.to_node in include:
            degree[link.from_node] = degree.get(link.from_node, 0) + 1
            degree[link.to_node] = degree.get(link.to_node, 0) + 1
            neighbors.setdefault(link.from_node, set()).add(link.to_node)
            neighbors.setdefault(link.to_node, set()).add(link.from_node)

    # Compute shared nodes from the full graph for coloring only.
    ref_path = _build_ref_nodes(g)
    shared_path = _build_shared_nodes(g)
    sample_counts = _sample_count_map(g)

    nodes = []
    for seg_name in include:
        seg = g.segments.get(seg_name)
        if seg is None:
            continue
        nodes.append({
            "id": seg_name,
            "length": seg.length,
            "on_selected_path": seg_name in selected,
            "is_shared": seg_name in shared_path,
            "on_reference_path": seg_name in ref_path,
            "sample_count": sample_counts.get(seg_name, 0),
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
            "is_shared": link.from_node in shared_path and link.to_node in shared_path,
        })

    original_nodes = len(nodes)
    original_edges = len(edges)
    truncated = original_nodes > WEB_GRAPH_MAX_NODES or original_edges > WEB_GRAPH_MAX_EDGES
    nodes = nodes[:WEB_GRAPH_MAX_NODES]
    # Filter edges to only those whose both endpoints are in the truncated node set.
    truncated_ids = set(n["id"] for n in nodes)
    edges = [e for e in edges if e["source"] in truncated_ids and e["target"] in truncated_ids]
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
        "stitched": ["results/merge/merged.gfa", "work/demo/merged.gfa"],
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
        "baseline": "AVAILABLE" if "baseline" in graphs else "NOT_AVAILABLE",
        "parallel_chunks": "IMPLEMENTED" if chunks else "NOT_AVAILABLE",
        "stitch": (
            "PLACEHOLDER_LINEAR_ONLY"
            if "stitched" in graphs and graph_meta.get("stitched", {}).get("paths", 0) <= 5
            and graph_meta.get("stitched", {}).get("edges", 0) < 1000
            else "IMPLEMENTED" if "stitched" in graphs
            else "NOT_IMPLEMENTED"
        ),
        "path_equivalence": "NOT_RUN",
        "variant_equivalence": "NOT_RUN",
    }

    has_real_baseline = pipeline_status["baseline"] == "AVAILABLE"
    has_real_stitch = pipeline_status["stitch"] not in ("NOT_IMPLEMENTED", "PLACEHOLDER_LINEAR_ONLY")
    if has_real_baseline and has_real_stitch:
        data_mode = "real"
    elif has_real_baseline:
        data_mode = "real_baseline_pending_stitch"
    else:
        data_mode = "synthetic"

    samples = _manifest_samples(samples_by_graph)
    target = _load_config_target()
    git_sha = _git_sha()
    generated_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    manifest = {
        "schema_version": "2",
        "run": {
            "run_id": "latest",
            "git_sha": git_sha,
            "generated_at": generated_at,
            "data_mode": data_mode,
            "scientific_status": (
                "Real HPRC baseline available. "
                "Parallel stitch pending chunk rebuild (mash-kmer=19)."
                if data_mode == "real_baseline_pending_stitch"
                else "Full real HPRC dataset."
                if data_mode == "real"
                else "Synthetic demo dataset."
            ),
        },
        "target": target,
        "graphs": graph_meta,
        "pipeline_status": pipeline_status,
        "samples": [{"sample": s["sample"],
                     "haplotypes": s["haplotypes"],
                     "hap_labels": {h: _hap_label(h) for h in s["haplotypes"]}}
                    for s in samples],
    }
    with open("results/web/manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    overview = {
        "schema_version": "2",
        "run": manifest["run"],
        "target": target,
        "chunks": chunks,
        "graphs": graph_meta,
        "pipeline_status": pipeline_status,
    }
    with open("results/web/overview.json", "w") as f:
        json.dump(overview, f, indent=2)

    # Path diagnostics
    diag_rows = []
    for label, path in graphs.items():
        g = GfaGraph.parse_file(path)
        for sample, haps in sorted(samples_by_graph[label].items()):
            for hap in sorted(haps):
                pn, steps = _build_selected_path(g, sample, hap)
                if not pn:
                    continue
                diag_rows.append({
                    "graph": label,
                    "sample": sample,
                    "haplotype": hap,
                    "path_name": pn,
                    "step_count": len(steps),
                    "length_bp": sum(g.segments[s].length
                                   for s, _ in steps if s in g.segments),
                    "gfa_source": path,
                })
    with open("results/web/path_diagnostics.json", "w") as f:
        json.dump(diag_rows, f, indent=2)
    for row in diag_rows:
        if row["step_count"] <= 3 and row["sample"] != "GRCh38":
            print(f"WARNING: {row['sample']}#{row['haplotype']} has only "
                  f"{row['step_count']} step(s) in {row['graph']} "
                  f"— inspect PGGB GFA / parser.", file=sys.stderr)

    print(f"build_web_dataset: {len(graphs)} graph(s), {len(samples)} sample(s), "
          f"{written} sample graph JSON(s), {len(chunks)} chunk(s)")
    print("  manifest  -> results/web/manifest.json")
    print("  overview  -> results/web/overview.json")
    print("  path diag -> results/web/path_diagnostics.json")


if __name__ == "__main__":
    main()
