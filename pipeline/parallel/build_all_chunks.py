"""build_all_chunks.py - Build per-chunk FASTA by independently aligning
each reference chunk against each HPRC assembly via minimap2.

NO linear scaling. NO loading entire genomes — uses samtools faidx on
the exact source_contig for each independently mapped chunk.
"""
import argparse, csv, hashlib, json, os, subprocess, sys, yaml

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


def file_sha256(path):
    if not os.path.exists(path): return ""
    with open(path, "rb") as f: return hashlib.sha256(f.read()).hexdigest()


def git_commit(path="."):
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                           text=True, cwd=path)
        return r.stdout.strip()[:12] if r.returncode == 0 else "no-git"
    except: return "no-git"


def compute_provenance(cid, cs, ce, cfg_p="config/pipeline.yaml",
                       map_p="results/preparation/sequence_mapping.tsv"):
    return hashlib.sha256(
        f"{git_commit()}|{file_sha256(cfg_p)}|{file_sha256(map_p)}|{cid}|{cs}|{ce}"
        .encode()).hexdigest()[:16]


def load_provenance(path="work/chunks/provenance.json"):
    if not os.path.exists(path): return {}
    with open(path) as f: return json.load(f)


def save_provenance(prov, path="work/chunks/provenance.json"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f: json.dump(prov, f, indent=2)
    return best if bc >= min_cov else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--execute", action="store_true")
    p.add_argument("--force", action="store_true",
                   help="Force rebuild even if provenance matches")
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

    # Build parent-locus index from sequence_mapping.tsv
    parent_locus = {}
    for r in mapping:
        n = r["assembly_name"]
        if n not in parent_locus:
            parent_locus[n] = {
                "sample": r["sample"], "haplotype": r["haplotype"],
                "haplotype_label": r["haplotype_label"],
                "source_contig": r.get("source_contig", ""),
                "source_start": int(r.get("source_start", 0) or 0),
                "source_end": int(r.get("source_end", 0) or 0),
                "strand": r.get("strand", "+"),
            }

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
    prov = load_provenance()
    stale_count = 0; fresh_count = 0
    print(f"Building {len(chunks)} alignment-projected chunk FASTA files...")

    for c in chunks:
        cid = c["chunk_id"]
        cs = int(c["reference_start"]); ce = int(c["reference_end"])
        op = f"work/chunks/{cid}.fa"
        prov_key = f"chunk_{cid}"
        expected_prov = compute_provenance(cid, cs, ce)

        # Provenance check — regenerate when stale
        if not args.force and os.path.exists(op) and os.path.getsize(op) > 0:
            if prov.get(prov_key) == expected_prov:
                print(f"  CACHED {cid} (provenance match)")
                fresh_count += 1
                continue
            else:
                print(f"  STALE {cid} (provenance mismatch, regenerating)")
                stale_count += 1
        else:
            fresh_count += 1

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

                # === PARENT-LOCUS CONSTRAINT ===
                pl = parent_locus.get(name, {})
                expected_ctg = pl.get("source_contig", "")
                if expected_ctg and ctg != expected_ctg:
                    print(f"    FATAL: {sm} ({hl}) chunk {cid} maps to {ctg} "
                          f"but parent 1Mb mapping is on {expected_ctg}. "
                          f"Chunk jumps to another contig/paralog!", file=sys.stderr)
                    sys.exit(1)
                actual_strand = res["strand"]
                expected_strand = pl.get("strand", "+")
                if actual_strand != expected_strand:
                    print(f"    WARNING: {sm} ({hl}) chunk {cid} strand {actual_strand} "
                          f"vs parent {expected_strand}. Using actual strand.")
                ps = pl.get("source_start", 0)
                pe = pl.get("source_end", 0)
                if ps and pe:
                    margin = max(ce - cs, 100000)
                    ts = res["source_start"]
                    te = res["source_end"]
                    if ts < ps - margin or te > pe + margin:
                        print(f"    WARNING: {sm} ({hl}) chunk {cid} source [{ts},{te}) "
                              f"outside parent [{ps},{pe}) + margin {margin}")
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
        prov[prov_key] = expected_prov

    save_provenance(prov)

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
    print(f"Done. {len(chunks)} chunks. Fresh: {fresh_count}, Stale: {stale_count}")

if __name__=="__main__": main()
