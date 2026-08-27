"""Compute graph statistics for baseline and merged GFAs.

Provides topology analysis (components, degrees, N50) used by
benchmark reports, the web dashboard, and merge validation.
"""
import json
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph

class _DSU:
    """Union-Find data structure for connected-component analysis."""
    def __init__(self, elements):
        self._parent = {e: e for e in elements}
        self._size = {e: 1 for e in elements}
    def find(self, x):
        while self._parent[x] != x:
            self._parent[x] = self._parent[self._parent[x]]
            x = self._parent[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb: return
        if self._size[ra] < self._size[rb]: ra, rb = rb, ra
        self._parent[rb] = ra
        self._size[ra] += self._size[rb]
    def components(self):
        roots = {}
        for e in self._parent:
            r = self.find(e)
            roots[r] = roots.get(r, 0) + 1
        return sorted(roots.values(), reverse=True)

def _n50(lengths):
    if not lengths: return 0
    sorted_l = sorted(lengths, reverse=True)
    total = sum(sorted_l)
    half = total / 2.0
    cum = 0
    for val in sorted_l:
        cum += val
        if cum >= half: return val
    return sorted_l[-1] if sorted_l else 0

LEX_TSV_COLUMNS = ["label","n50","nodes","edges","paths","walks",
    "component_count","seqwish_ratio","sequences","status",
    "wall_s","peak_mib","threads"]

def compute_stats(graph, label=""):
    nodes = graph.node_count()
    edges = graph.edge_count()
    seg_names = list(graph.segments.keys())
    dsu = _DSU(seg_names)
    for link in graph.links:
        if link.from_node in graph.segments and link.to_node in graph.segments:
            dsu.union(link.from_node, link.to_node)
    comps = dsu.components()
    comp_count = len(comps)
    comp_largest = comps[0] if comps else 0
    clf = round(comp_largest / max(nodes, 1), 4) if nodes else 0.0
    in_deg, out_deg = {}, {}
    seen, dangling = set(), 0
    for link in graph.links:
        out_deg[link.from_node] = out_deg.get(link.from_node, 0) + 1
        in_deg[link.to_node] = in_deg.get(link.to_node, 0) + 1
        seen.add((link.from_node, link.from_orient, link.to_node, link.to_orient))
        if link.from_node not in graph.segments or link.to_node not in graph.segments:
            dangling += 1
    branching = sum(1 for n in seg_names
                    if max(out_deg.get(n,0), in_deg.get(n,0)) > 1)
    max_out = max(out_deg.values()) if out_deg else 0
    max_in = max(in_deg.values()) if in_deg else 0
    tips = sum(1 for n in seg_names
               if (in_deg.get(n,0)==0 or out_deg.get(n,0)==0))
    uniq = len(seen)
    paths = graph.path_count()
    walks = graph.walk_count()
    n_paths = paths + walks
    tss = 0; pp = {}
    for p in graph.paths.values():
        st = len(p.segment_names) if p.segment_names else 0
        tss += st; pp[p.path_name] = st
    for w in graph.walks:
        st = len(w.path) if w.path else 0
        tss += st
        pan = w.sample + "#" + w.haplotype + "#" + w.contig
        pp[pan] = pp.get(pan, 0) + st
    tb = graph.total_sequence_bp()
    sl = [s.length for s in graph.segments.values()]
    enr = round(edges / max(nodes,1), 4)
    bf = round(branching / max(nodes,1), 4) if nodes else 0.0
    return {
        "label": label, "nodes": nodes, "edges": edges,
        "node_lengths": {"total_bp": tb, "n50": _n50(sl)},
        "components": {"count": comp_count, "largest": comp_largest,
                        "largest_fraction": clf},
        "degrees": {"branching_nodes": branching, "max_out_degree": max_out,
                     "max_in_degree": max_in, "tip_nodes": tips,
                     "unique_edges": uniq, "dangling_links": dangling},
        "paths": paths, "walks": walks,
        "path_steps": {"count": n_paths, "total_steps": tss,
                         "per_path": dict(sorted(pp.items()))},
        "complexity": {"edge_node_ratio": enr, "branching_fraction": bf},
        "total_bp": tb,
        "samples": sorted(graph.get_sample_names()),
    }

def compare(a, b, tolerance_pct=5.0):
    metrics = ["nodes","edges","paths","walks","component_count","paths_plus_walks"]
    checks = []; n_fail = 0
    for m in metrics:
        if m == "component_count":
            va = a.get("components",{}).get("count",0)
            vb = b.get("components",{}).get("count",0)
        elif m == "paths_plus_walks":
            va = a.get("paths",0) + a.get("walks",0)
            vb = b.get("paths",0) + b.get("walks",0)
        else:
            va, vb = a, b
            for k in m.split("."):
                va = va.get(k, 0) if isinstance(va, dict) else 0
                vb = vb.get(k, 0) if isinstance(vb, dict) else 0
        if not isinstance(va, (int,float)): va = 0
        if not isinstance(vb, (int,float)): vb = 0
        if va == 0 and vb == 0:
            dp = 0.0; strict = False; status = "PASS"
        elif va == 0:
            dp = 100.0; strict = True; status = "FAIL"
        else:
            dp = round(abs(va-vb)/max(va,1)*100, 2)
            strict = dp > tolerance_pct
            status = "PASS" if not strict else "FAIL"
        checks.append({"metric":m,"baseline":va,"merged":vb,
            "delta_pct":dp,"strict":strict,"status":status})
        if strict: n_fail += 1
    ca = a.get("components",{}).get("count",0)
    cb = b.get("components",{}).get("count",0)
    eq = ca == cb
    v = "EQUIVALENT" if eq and n_fail == 0 else "DIVERGENT"
    return {"verdict":v,"n_fail":n_fail,"checks":checks}

def web_metrics_block(all_stats):
    block = {}
    for label, s in all_stats.items():
        if label == "comparison": continue
        block[label] = {"nodes":s.get("nodes",0),
            "edges":s.get("edges",0),
            "paths":s.get("paths",0)+s.get("walks",0)}
    return block

def _results_dir_from_config(config_path):
    try:
        import yaml
        with open(config_path) as f: cfg = yaml.safe_load(f)
        return cfg.get("output",{}).get("results_dir","results")
    except Exception: return "results"

def append_lex_tsv_row(stats_dict, tsv_path, **kwargs):
    s = stats_dict
    row = [str(s.get("label",""))]
    row.append(str(s.get("node_lengths",{}).get("n50",0)))
    for k in ("nodes","edges","paths","walks"):
        row.append(str(s.get(k,0)))
    row.append(str(s.get("components",{}).get("count",0)))
    er = round(s.get("complexity",{}).get("edge_node_ratio",0),4)
    row.append(str(er))
    row.append(str(s.get("paths",0) + s.get("walks",0)))
    row.append("OK")
    for k in ("wall_s","peak_mib","threads"):
        row.append(str(kwargs.get(k,"")))
    line = "\t".join(row) + "\n"
    if not os.path.exists(tsv_path):
        with open(tsv_path,"w") as fh:
            fh.write("\t".join(LEX_TSV_COLUMNS)+"\n")
    with open(tsv_path,"a") as fh: fh.write(line)

def stats(g, label=""):
    return dict(label=label, nodes=g.node_count(), edges=g.edge_count(),
                paths=g.path_count(), walks=g.walk_count(),
                total_bp=g.total_sequence_bp(),
                samples=sorted(g.get_sample_names()))

def main():
    rd = "results"
    targets = [(f"{rd}/baseline/baseline.gfa","baseline"),
               (f"{rd}/merge/merged.gfa","merged")]
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

if __name__ == "__main__": main()
