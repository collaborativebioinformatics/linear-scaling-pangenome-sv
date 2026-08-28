"""Synthetic end-to-end demo: sequences, GFAs, merge, JSON."""
import json
import os
import sys
import yaml
import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipeline.merge.gfa import GfaGraph, Header, Segment, Link, Path
from pipeline.merge.merge_graphs import (
    diagnostic_disjoint_union, overlap_aware_stitch, _write_boundary_report,
)
from pipeline.merge.validate_merge import validate
from pipeline.benchmark.graph_stats import compute_stats, compare
from pipeline.parallel.make_chunks import create_chunks, write_manifest

SEG = 20


def _chain(g, seq, label):
    prev, nodes = None, []
    for i in range(0, len(seq), SEG):
        e = min(i + SEG, len(seq))
        nn = f"{label}_s{i // SEG}"
        g.segments[nn] = Segment(nn, seq[i:e])
        if prev:
            g.links.append(Link(prev, "+", nn, "+"))
        prev = nn
        nodes.append(nn)
    return nodes


def build_sequences(cfg):
    rl = cfg["demo"]["reference_length"]
    bp = "GGGGGAAAAACCCCCTTTTT"
    ref = (bp * (rl // 20 + 1))[:rl]
    rp = cfg["demo"]["output"]["reference_fasta"]
    os.makedirs(os.path.dirname(rp), exist_ok=True)
    with open(rp, "w") as f:
        f.write(">GRCh38#0#chr21\n")
        for i in range(0, len(ref), 80):
            f.write(ref[i:i + 80] + "\n")
    configs = [
        ("H1", "HG00673", "1", 400, "ins", "AGCT" * 13),
        ("H2", "HG00673", "2", 800, "del", 30),
        ("H3", "HG00733", "1", 1200, "snp", None),
        ("H4", "HG00733", "2", 1600, "complex", "CCCCAAAA"),
    ]
    haps = []
    hp = cfg["demo"]["output"]["haplotypes_fasta"]
    with open(hp, "w") as f:
        for _n, sm, hap, pos, vt, vd in configs:
            seq = ref
            if vt == "ins" and vd:
                seq = ref[:pos] + vd + ref[pos:]
            elif vt == "del":
                seq = ref[:pos] + ref[pos + 30:]
            elif vt == "snp":
                sl = list(seq)
                for i in range(pos, min(pos + 5, len(sl))):
                    sl[i] = "A" if sl[i] in "GC" else "T"
                seq = "".join(sl)
            elif vt == "complex" and vd:
                seq = ref[:pos] + vd + ref[pos + len(vd):]
            f.write(f">{sm}#{hap}#chr21 {_n}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
            haps.append({"sample": sm, "haplotype": hap, "sequence": seq})
    return ref, haps


def baseline(ref, haps, cfg):
    g = GfaGraph()
    g.headers.append(Header("1.1"))
    for seq, sm, hap in [(ref, "GRCh38", "0")] + \
                         [(h["sequence"], h["sample"], h["haplotype"]) for h in haps]:
        nodes = _chain(g, seq, f"bl_{sm}_{hap}")
        if nodes:
            g.paths[f"{sm}#{hap}#chr21"] = Path(
                f"{sm}#{hap}#chr21", [n + "+" for n in nodes])
    op = cfg["demo"]["output"]["baseline_gfa"]
    os.makedirs(os.path.dirname(op), exist_ok=True)
    g.write_gfa(op)
    return g


def chunks(ref, haps, cfg):
    cc = cfg["demo"]["chunk"]
    cd = cfg["demo"]["output"]["chunk_dir"]
    os.makedirs(cd, exist_ok=True)
    os.makedirs("work/demo", exist_ok=True)
    chs = create_chunks("chr21", 0, len(ref), cc["size"], cc["overlap"])
    write_manifest(chs, "work/demo/chunk_manifest.tsv")
    result = []
    for c in chs:
        cid = c["chunk_id"]
        st, en = c["reference_start"], c["reference_end"]
        cg = GfaGraph()
        cg.headers.append(Header("1.1"))
        rn = _chain(cg, ref[st:en], f"{cid}_ref")
        if rn:
            rpn = f"GRCh38#0#chr21:{st}-{en}"
            cg.paths[rpn] = Path(rpn, [n + "+" for n in rn])
        last = c is chs[-1]
        for h in haps:
            hs, he = st, min(en, len(h["sequence"]))
            if last:
                he = len(h["sequence"])
            if he <= hs:
                continue
            hn = _chain(cg, h["sequence"][hs:he],
                        f"{cid}_{h['sample']}_{h['haplotype']}")
            if hn:
                pn = f"{h['sample']}#{h['haplotype']}#chr21:{hs}-{he}"
                cg.paths[pn] = Path(pn, [n + "+" for n in hn])
        cg.write_gfa(os.path.join(cd, f"{cid}.gfa"))
        result.append((cid, cg))
    return result, chs


def main():
    cfg = yaml.safe_load(open("config/demo.yaml"))
    pipe = {}
    if os.path.exists("config/pipeline.yaml"):
        pipe = yaml.safe_load(open("config/pipeline.yaml")) or {}
    strategy = pipe.get("merge", {}).get("strategy", "overlap_aware")
    print("Demo: building synthetic data...")
    print(f"Strategy: {strategy}")
    ref, haps = build_sequences(cfg)
    bg = baseline(ref, haps, cfg)
    cgs, chs = chunks(ref, haps, cfg)
    rows = {c["chunk_id"]: c for c in chs}
    obp = cfg["demo"]["chunk"]["overlap"]

    if strategy == "disjoint_union":
        merged = diagnostic_disjoint_union(cgs)
        br = []
    else:
        merged, br = overlap_aware_stitch(
            cgs, overlap_bp=obp, chunk_rows=rows)

    mp = cfg["demo"]["output"]["merged_gfa"]
    os.makedirs(os.path.dirname(mp), exist_ok=True)
    merged.write_gfa(mp)
    os.makedirs("results/merge", exist_ok=True)
    merged.write_gfa("results/merge/merged.gfa")
    if br:
        _write_boundary_report(br, "work/demo/boundary_report.tsv")
        _write_boundary_report(br, "results/merge/boundary_report.tsv")

    issues = validate(cfg["demo"]["output"]["baseline_gfa"], mp)
    cmp = compare(compute_stats(bg, "baseline"), compute_stats(merged, "merged"))
    comps = compute_stats(merged, "merged")["components"]["count"]
    stitch_ok = bool(br) and all(b.get("status") == "PASS" for b in br)
    print(f"Baseline: {bg.node_count()}n {bg.edge_count()}e")
    print(f"Merged: {merged.node_count()}n {merged.edge_count()}e "
          f"{merged.path_count()} paths, {comps} components")
    print(f"Verdict: {cmp['verdict']}")
    if br:
        failed = [b for b in br if b["status"] == "FAIL"]
        print(f"Boundaries: {len(br) - len(failed)}/{len(br)} PASS")
    errors = [i for i in issues if i["severity"] == "ERROR"]
    if errors:
        for e in errors[:5]:
            print(f"  [ERROR] {e['message']}")

    data = {
        "data_mode": "synthetic",
        "run": {"run_id": "demo", "pipeline_version": "0.1.0", "mode": "demo",
                "data_mode": "synthetic", "strategy": strategy,
                "timestamp": datetime.datetime.now().isoformat()},
        "target": {"reference": "GRCh38", "chromosome": "chr21",
                   "start": 0, "end": cfg["demo"]["reference_length"]},
        "samples": ["HG00673", "HG00733"],
        "metrics": {
            "baseline": {"nodes": bg.node_count(), "edges": bg.edge_count()},
            "merged": {"nodes": merged.node_count(), "edges": merged.edge_count(),
                       "paths": merged.path_count(), "components": comps},
        },
        "stitch": {"status": "PASS" if stitch_ok else ("FAIL" if br else "OK"),
                   "strategy": strategy},
        "equivalence": {"verdict": cmp["verdict"], "n_fail": cmp["n_fail"],
                        "checks": cmp["checks"]},
        "boundaries": br, "bubbles": [], "graphWindow": {},
        "validation_errors": len(errors),
    }
    jp = cfg["demo"]["output"]["json_output"]
    os.makedirs(os.path.dirname(jp), exist_ok=True)
    json.dump(data, open(jp, "w"), indent=2)
    print(f"JSON: {jp}")
    print("Run: cd web && npm run dev")


if __name__ == "__main__":
    main()