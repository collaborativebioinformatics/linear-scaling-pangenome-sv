"""build_all_chunks.py - Build per-chunk FASTA by independently aligning
each reference chunk against each HPRC assembly via minimap2.

NO linear scaling. NO loading entire genomes — uses samtools faidx on
the exact source_contig for each independently mapped chunk.
"""
import argparse, csv, os, subprocess, sys, yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.prepare.faidx_utils import faidx_extract, ensure_faidx, _revcomp

def load_manifest(path="work/chunks/chunk_manifest.tsv"):
    if not os.path.exists(path): return []
    with open(path) as f: return list(csv.DictReader(f, delimiter="\t"))

def load_mapping(path="results/preparation/sequence_mapping.tsv"):
    if not os.path.exists(path): return []
    return list(csv.DictReader(open(path), delimiter="\t"))

def find_assembly_path(name):
    for d in ["work/downloads","/data/hprc"]:
        for ext in [".fa.gz",".fa",".fasta",".fna"]:
            p = os.path.join(d, name+ext)
            if os.path.exists(p): return p
    return None

def extract_chunk_fasta(ref_seq, cs, ce, out):
    c = ref_seq[cs:ce]
    with open(out,"w") as f:
        f.write(f">GRCh38#0#chr21:{cs}-{ce}\n")
        for i in range(0,len(c),80): f.write(c[i:i+80]+"\n")

def map_chunk(ap, qry, min_mapq, min_cov):
    if subprocess.run(["which","minimap2"],capture_output=True).returncode != 0: return None
    try:
        r = subprocess.run(["minimap2","-x","asm5","-t","4","--eqx","-c",str(ap),qry],
                           capture_output=True,text=True,timeout=600)
    except: return None
    best, bc = None, 0
    for line in r.stdout.strip().split("\n"):
        if not line: continue
        p = line.split("\t")
        if len(p) < 12: continue
        qlen=int(p[1]); qs=int(p[2]); qe=int(p[3])
        st=p[4]; tn=p[5]; ts=int(p[7]); te=int(p[8])
        mat=int(p[9]); bl=int(p[10]); mpq=int(p[11])
        qc=(qe-qs)/max(qlen,1)
        if mpq >= min_mapq and qc > bc:
            bc=qc
            best={"source_contig":tn,"source_start":ts,"source_end":te,"strand":st,
                  "mapping_quality":mpq,"query_coverage":qc,
                  "identity":mat/max(bl,1) if bl>0 else 0,
                  "method":"minimap2_asm5_chunk"}
    return best if bc >= min_cov else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--execute", action="store_true")
    args = p.parse_args()
    chunks = load_manifest()
    if not chunks: print("No manifest."); return
    mapping = load_mapping()
    if not mapping: print("No mapping."); return

    # Read mapping thresholds from config
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    mc = cfg.get("mapping", {})
    min_mapq = mc.get("min_mapq", 20)
    min_cov = mc.get("min_query_coverage", 0.90)

    # Reference — single contig, read into memory once (~47 Mb)
    ref = "work/reference/GRCh38_chr21.fa"
    if not os.path.exists(ref): print("Ref not found."); return
    with open(ref) as f:
        ref_seq = "".join(line.strip() for line in f if not line.startswith(">"))

    # Pre-find HPRC assembly paths (do NOT load sequences into RAM)
    hap_paths = {}
    for r in mapping:
        n = r["assembly_name"]
        if n in hap_paths: continue
        ap = find_assembly_path(n)
        if ap:
            hap_paths[n] = ap
            ensure_faidx(ap)  # ensure indexed
        else:
            print(f"FATAL: {n} not found"); sys.exit(1)

    os.makedirs("work/chunks", exist_ok=True)
    os.makedirs("results/preparation", exist_ok=True)
    cm_rows = []
    print(f"Building {len(chunks)} alignment-projected chunk FASTA files...")

    for c in chunks:
        cid = c["chunk_id"]
        cs = int(c["reference_start"]); ce = int(c["reference_end"])
        op = f"work/chunks/{cid}.fa"
        if os.path.exists(op) and os.path.getsize(op) > 0:
            print(f"  EXISTS {cid}"); continue

        qf = f"work/preparation/{cid}_query.fa"
        extract_chunk_fasta(ref_seq, cs, ce, qf)
        sc = 0
        with open(op,"w") as fout:
            rc = ref_seq[cs:ce]
            fout.write(">GRCh38#0#chr21\n")
            for i in range(0,len(rc),80): fout.write(rc[i:i+80]+"\n")
            sc += 1
            for r in mapping:
                sm = r["sample"]; nh = r["haplotype"]
                hl = r["haplotype_label"]; name = r["assembly_name"]
                ap = hap_paths[name]
                res = map_chunk(ap, qf, min_mapq, min_cov)
                if not res:
                    print(f"    FATAL: {sm} ({hl}) chunk {cid} unmapped", file=sys.stderr)
                    sys.exit(1)
                ctg = res["source_contig"]
                ss = res["source_start"]; se = res["source_end"]
                strand = res["strand"]
                # Extract from exact source_contig via samtools faidx
                # NEVER concatenate all contigs — PAF coords are per-contig
                seq = faidx_extract(ap, ctg, ss, se, strand)
                fout.write(f">{sm}#{nh}#chr21\n")
                for i in range(0,len(seq),80): fout.write(seq[i:i+80]+"\n")
                sc += 1
                cm_rows.append({"chunk_id":cid,"sample":sm,"haplotype":nh,
                    "haplotype_label":hl,"assembly_name":name,
                    "ref_start":cs,"ref_end":ce,
                    "source_contig":ctg,"source_start":ss,"source_end":se,
                    "strand":strand,"identity":f"{res.get('identity',0):.4f}",
                    "query_coverage":f"{res.get('query_coverage',0):.4f}",
                    "mapping_quality":res["mapping_quality"],
                    "mapping_method":res["method"]})
        if os.path.exists(qf): os.remove(qf)
        print(f"  {cid}: {cs}-{ce} ({sc} seqs) -> {op}")

    if cm_rows:
        cm_path = "results/preparation/chunk_mapping.tsv"
        fn = ["chunk_id","sample","haplotype","haplotype_label","assembly_name",
              "ref_start","ref_end","source_contig","source_start","source_end",
              "strand","identity","query_coverage","mapping_quality","mapping_method"]
        with open(cm_path,"w",newline="") as f:
            w = csv.DictWriter(f,fieldnames=fn,delimiter="\t")
            w.writeheader(); w.writerows(cm_rows)
        print(f"  Chunk mapping: {cm_path} ({len(cm_rows)} records)")

    if args.execute:
        print("\nBuilding chunk graphs via Docker...")
        for c in chunks:
            cid = c["chunk_id"]
            if not os.path.exists(f"work/chunks/{cid}.gfa"):
                subprocess.run(["bash","pipeline/parallel/build_chunk.sh",cid], check=False)
    print(f"Done. {len(chunks)} chunks.")

if __name__=="__main__": main()
