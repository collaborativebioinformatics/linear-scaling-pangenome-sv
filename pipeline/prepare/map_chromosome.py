"""
map_chromosome.py - Map the GRCh38 smoke interval to each HPRC assembly
via minimap2, recording source_start/source_end/strand/query_coverage.

Fails hard if any interval fails MAPQ or query coverage thresholds.
"""
import csv
import gzip
import os
import subprocess
import sys
import yaml


def _revcomp(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G",
            "a": "t", "t": "a", "g": "c", "c": "g",
            "N": "N", "n": "n"}
    return "".join(comp.get(c, c) for c in reversed(seq))


def _open_read(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def extract_interval(fa, start, end, out, header="query"):
    """Extract a 0-based half-open interval from a FASTA file."""
    seq = []
    with _open_read(fa) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    s = "".join(seq)[start:end]
    with open(out, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(s), 80):
            f.write(s[i:i + 80] + "\n")


def map_interval(ap, qry, min_mapq, min_cov):
    """Align the query interval against the assembly via minimap2 asm5."""
    if subprocess.run(["which", "minimap2"],
                       capture_output=True).returncode != 0:
        return None
    try:
        r = subprocess.run(
            ["minimap2", "-x", "asm5", "-t", "4", "--eqx", "-c",
             str(ap), qry],
            capture_output=True, text=True, timeout=600)
    except Exception:
        return None
    best, bc = None, 0
    for line in r.stdout.strip().split("\n"):
        if not line:
            continue
        p = line.split("\t")
        if len(p) < 12:
            continue
        qlen = int(p[1])
        qs = int(p[2])
        qe = int(p[3])
        st = p[4]
        tn = p[5]
        ts = int(p[7])
        te = int(p[8])
        mpq = int(p[11])
        qc = (qe - qs) / max(qlen, 1)
        if mpq >= min_mapq and qc > bc:
            bc = qc
            best = {
                "source_contig": tn,
                "source_start": ts,
                "source_end": te,
                "strand": st,
                "mapping_quality": mpq,
                "query_coverage": qc,
                "method": "minimap2_asm5",
            }
    return best if bc >= min_cov else None


def main():
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    tgt = cfg["target"]
    mapping_cfg = cfg.get("mapping", {})
    min_mapq = mapping_cfg.get("min_mapq", 20)
    min_cov = mapping_cfg.get("min_query_coverage", 0.90)

    chrom, rs, re = tgt["chromosome"], tgt["start"], tgt["end"]
    print(f"=== Interval Mapping: {chrom}:{rs}-{re} "
          f"({re - rs} bp, 0-based half-open) ===")
    print(f"  Thresholds: MAPQ>={min_mapq}, query_coverage>={min_cov}")

    mp = "work/manifests/hprc_selected.csv"
    if not os.path.exists(mp):
        print("No manifest."); sys.exit(1)
    samples = list(csv.DictReader(open(mp)))

    os.makedirs("results/preparation", exist_ok=True)
    ref = "work/reference/GRCh38_chr21.fa"
    if not os.path.exists(ref):
        print("Ref not found."); sys.exit(1)

    qf = "work/preparation/smoke_interval_query.fa"
    extract_interval(ref, rs, re, qf,
                     f"GRCh38#0#{chrom}:{rs}-{re}")

    rows, failed = [], []
    for s in samples:
        sm = s.get("sample_id", "?")
        nh = s.get("haplotype", "?")
        hl = s.get("haplotype_label", "?")
        name = s.get("assembly_name", "?")
        print(f"  {sm} hap={nh} ({hl}): {name}")

        ap = None
        for d in ["work/downloads", "/data/hprc"]:
            for ext in [".fa.gz", ".fa", ".fasta", ".fna"]:
                p = os.path.join(d, name + ext)
                if os.path.exists(p):
                    ap = p
                    break
            if ap:
                break
        if not ap:
            failed.append(f"{sm} ({hl}): file not found")
            continue

        res = map_interval(ap, qf, min_mapq, min_cov)
        if not res:
            failed.append(f"{sm} ({hl}): interval "
                          f"{chrom}:{rs}-{re} unmapped or below threshold")
            continue

        print(f"    -> {res['source_contig']}:{res['source_start']}-"
              f"{res['source_end']} ({res['strand']}, "
              f"mapq={res['mapping_quality']}, "
              f"cov={res['query_coverage']:.3f})")

        rows.append({
            "sample": sm, "haplotype": nh,
            "haplotype_label": hl, "assembly_name": name,
            "reference_chromosome": chrom,
            "reference_start": rs, "reference_end": re,
            "source_contig": res["source_contig"],
            "source_start": res["source_start"],
            "source_end": res["source_end"],
            "strand": res["strand"],
            "mapping_method": res["method"],
            "mapping_quality": res["mapping_quality"],
            "query_coverage": f"{res['query_coverage']:.4f}",
            "confidence": "high" if res["mapping_quality"] >= 40
                          else "moderate",
            "status": "mapped",
        })

    mo = "results/preparation/sequence_mapping.tsv"
    fn = ["sample", "haplotype", "haplotype_label", "assembly_name",
          "reference_chromosome", "reference_start", "reference_end",
          "source_contig", "source_start", "source_end", "strand",
          "mapping_method", "mapping_quality", "query_coverage",
          "confidence", "status"]
    with open(mo, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fn, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\n{mo}: {len(rows)} records")
    if failed:
        print("\nFATAL:", file=sys.stderr)
        for f in failed:
            print(f"  {f}", file=sys.stderr)
        sys.exit(1)
    print("All 4 haplotypes interval-mapped. Proceeding.")


if __name__ == "__main__":
    main()