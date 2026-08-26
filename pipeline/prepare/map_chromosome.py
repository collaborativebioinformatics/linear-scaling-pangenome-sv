"""
map_chromosome.py — Map de novo assembly contigs to reference chromosomes.

SAFETY: Never assume GRCh38 coordinates match de novo assemblies.
Strategy: PanSN names, header scan, minimap2 fallback.
"""
import csv
import gzip
import os
import subprocess
import sys
import yaml


def _is_chr21(name):
    u = name.upper().strip()
    return u in ("CHR21", "21") or "CHR21" in u


def scan_contig_names(fasta_path):
    """Find chr21 contig via PanSN or direct header match."""
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    candidates = []
    with opener(fasta_path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            h = line[1:].strip().split()[0]
            parts = h.split("#")
            if len(parts) >= 3 and _is_chr21(parts[-1]):
                candidates.append({
                    "contig": h, "method": "PanSN",
                    "confidence": "high", "strand": "+"
                })
            elif _is_chr21(h):
                candidates.append({
                    "contig": h, "method": "header_match",
                    "confidence": "high", "strand": "+"
                })
    return candidates


def map_via_minimap2(assembly_path, ref_path):
    """Minimap2 fallback for contig-to-chromosome mapping."""
    if subprocess.run(["which", "minimap2"],
                       capture_output=True).returncode != 0:
        return None
    if not os.path.exists(ref_path):
        return None
    try:
        r = subprocess.run(
            ["minimap2", "-x", "asm5", "-t", "4", ref_path, str(assembly_path)],
            capture_output=True, text=True, timeout=300)
        best, best_cov = None, 0
        for line in r.stdout.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t")
            if len(p) < 10:
                continue
            cov = int(p[9]) / max(int(p[1]), 1)
            if cov > best_cov:
                best_cov = cov
                best = {
                    "contig": p[0], "strand": p[4],
                    "method": "minimap2",
                    "confidence": "high" if cov > 0.5 else "moderate"
                }
        return best
    except Exception:
        return None


def main():
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    chrom = cfg["target"]["chromosome"]
    print("=== Chromosome Mapping ===")
    mp = "work/manifests/hprc_selected.csv"
    if not os.path.exists(mp):
        print("No manifest."); sys.exit(1)
    os.makedirs("results/preparation", exist_ok=True)
    samples = list(csv.DictReader(open(mp)))
    rows = []
    for s in samples:
        sm = s.get("sample_id", "?")
        hap = s.get("haplotype", "?")
        name = s.get("assembly_name", "?")
        print(f"  {sm} ({hap}): {name}")
        ap = None
        for d in ["work/downloads", "/data/hprc"]:
            for ext in ["", ".fa", ".fasta", ".fna", ".fa.gz"]:
                p = os.path.join(d, name + ext)
                if os.path.exists(p):
                    ap = p; break
            if ap:
                break
        if not ap:
            print(f"    NOT FOUND")
            rows.append({"sample": sm, "haplotype": hap,
                         "reference_chromosome": chrom,
                         "source_contig": "NOT_FOUND", "strand": "?",
                         "mapping_method": "N/A", "confidence": "N/A",
                         "status": "unresolved"})
            continue
        cands = scan_contig_names(ap)
        if cands:
            b = cands[0]
            print(f"    -> {b['contig']} ({b['method']})")
            rows.append({"sample": sm, "haplotype": hap,
                         "reference_chromosome": chrom,
                         "source_contig": b["contig"], "strand": b["strand"],
                         "mapping_method": b["method"],
                         "confidence": b["confidence"],
                         "status": "mapped"})
            continue
        mm = map_via_minimap2(ap, "work/reference/GRCh38_chr21.fa")
        if mm:
            print(f"    -> {mm['contig']} ({mm['method']})")
            rows.append({"sample": sm, "haplotype": hap,
                         "reference_chromosome": chrom,
                         "source_contig": mm["contig"], "strand": mm["strand"],
                         "mapping_method": mm["method"],
                         "confidence": mm["confidence"],
                         "status": "mapped"})
        else:
            print(f"    FAILED")
            rows.append({"sample": sm, "haplotype": hap,
                         "reference_chromosome": chrom,
                         "source_contig": "UNRESOLVED", "strand": "?",
                         "mapping_method": "failed", "confidence": "none",
                         "status": "unresolved"})
    mo = "results/preparation/sequence_mapping.tsv"
    with open(mo, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "sample", "haplotype", "reference_chromosome",
            "source_contig", "strand", "mapping_method",
            "confidence", "status"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    mapped = sum(1 for r in rows if r["status"] == "mapped")
    print(f"\n{mo}: {len(rows)} total, {mapped} mapped")
    unresolved = [r for r in rows
                  if r["status"] != "mapped"
                  and r["source_contig"] != "NOT_FOUND"]
    if unresolved:
        print("WARNING: Some chr21 mappings unresolved. "
              "Multi-FASTA will include full assemblies.")


if __name__ == "__main__":
    main()