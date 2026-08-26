"""Synthetic end-to-end demo: sequences, GFAs, merge, JSON."""
import json
import os
import sys
import yaml
import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipeline.merge.gfa import GfaGraph, Header, Segment, Link, Path
from pipeline.merge.merge_graphs import diagnostic_disjoint_union
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
            cg.paths[f"GRCh38#0#chr21_{cid}"] = Path(
                f"GRCh38#0#chr21_{cid}", [n + "+" for n in rn])
        for h in haps:
            hn = _chain(cg, h["sequence"][st:en],
                        f"{cid}_{h['sample']}_{h['haplotype']}")
            if hn:
                pn = f"{h['sample']}#{h['haplotype']}#chr21_{cid}"
                cg.paths[pn] = Path(pn, [n + "+" for n in hn])
        cg.write_gfa(os.path.join(cd, f"{cid}.gfa"))
        result.append((cid, cg))
    return result


def main():
    cfg = yaml.safe_load(open("config/demo.yaml"))
    print("Demo: building synthetic data...")
    ref, haps = build_sequences(cfg)
    bg = baseline(ref, haps, cfg)
    cgs = chunks(ref, haps, cfg)
    merged = diagnostic_disjoint_union(cgs)
    mp = cfg["demo"]["output"]["merged_gfa"]
    os.makedirs(os.path.dirname(mp), exist_ok=True)
    merged.write_gfa(mp)
    data = {
        "data_mode": "synthetic",
        "run": {"run_id": "demo", "pipeline_version": "0.1.0", "mode": "demo",
                "data_mode": "synthetic",
                "timestamp": datetime.datetime.now().isoformat()},
        "target": {"reference": "GRCh38", "chromosome": "chr21",
                   "start": 0, "end": cfg["demo"]["reference_length"]},
        "samples": ["HG00673", "HG00733"],
        "metrics": {
            "baseline": {"nodes": bg.node_count(), "edges": bg.edge_count()},
            "merged": {"nodes": merged.node_count(), "edges": merged.edge_count()},
        },
        "boundaries": [], "bubbles": [], "graphWindow": {},
    }
    jp = cfg["demo"]["output"]["json_output"]
    os.makedirs(os.path.dirname(jp), exist_ok=True)
    json.dump(data, open(jp, "w"), indent=2)
    print(f"Baseline: {bg.node_count()}n {bg.edge_count()}e")
    print(f"Merged: {merged.node_count()}n {merged.edge_count()}e")
    print(f"JSON: {jp}")
    print("Run: cd web && npm run dev")


if __name__ == "__main__":
    main()