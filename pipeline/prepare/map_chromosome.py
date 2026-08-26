"""
map_chromosome.py — Map de novo assembly contigs to reference chromosomes.

SAFETY: Never assume GRCh38 coordinates match de novo assemblies.
Strategy: PanSN names, header scan, minimap2 fallback.

FAILS HARD if any required haplotype cannot be mapped.
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
    failed = []
    for s in samples:
        sm = s.get("sample_id", "?")
        numeric_hap = s.get("haplotype", "?")
        hap_label = s.get("haplotype_label", "?")
        name = s.get("assembly_name", "?")
        print(f"  {sm} haplotype={numeric_hap} ({hap_label}): {name}")

        ap = None
        for d in ["work/downloads", "/data/hprc"]:
            for ext in [".fa.gz", ".fa", ".fasta", ".fna"]:
                p = os.path.join(d, name + ext)
                if os.path.exists(p):
                    ap = p; break
            if ap:
                break
        if not ap:
            print(f"    NOT FOUND — assembly file missing")
            failed.append(f"{sm} ({hap_label}): file {name} not found")
            continue

        cands = scan_contig_names(ap)
        if cands:
            b = cands[0]
            print(f"    -> chr21: {b['contig']} ({b['method']})")
            rows.append({"sample": sm, "haplotype": numeric_hap,
                         "haplotype_label": hap_label,
                         "assembly_name": name,
                         "reference_chromosome": chrom,
                         "source_contig": b["contig"], "strand": b["strand"],
                         "mapping_method": b["method"],
                         "confidence": b["confidence"],
                         "status": "mapped"})
            continue

        mm = map_via_minimap2(ap, "work/reference/GRCh38_chr21.fa")
        if mm:
            print(f"    -> chr21: {mm['contig']} ({mm['method']})")
            rows.append({"sample": sm, "haplotype": numeric_hap,
                         "haplotype_label": hap_label,
                         "assembly_name": name,
                         "reference_chromosome": chrom,
                         "source_contig": mm["contig"], "strand": mm["strand"],
                         "mapping_method": mm["method"],
                         "confidence": mm["confidence"],
                         "status": "mapped"})
        else:
            print(f"    FAILED — cannot determine chr21 contig")
            failed.append(f"{sm} ({hap_label}): chr21 mapping unresolved")

    # Write mapping report
    mo = "results/preparation/sequence_mapping.tsv"
    with open(mo, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "sample", "haplotype", "haplotype_label", "assembly_name",
            "reference_chromosome", "source_contig", "strand",
            "mapping_method", "confidence", "status"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\n{mo}: {len(rows)} records written")

    # FAIL HARD if any haplotype could not be mapped
    if failed:
        print("\nFATAL: The following required haplotypes could not be mapped to chr21:",
              file=sys.stderr)
        for f in failed:
            print(f"  {f}", file=sys.stderr)
        print("\nPipeline cannot proceed without valid chr21 mapping for all 4 haplotypes.",
              file=sys.stderr)
        sys.exit(1)

    print("All 4 haplotypes successfully mapped to chr21. Proceeding.")


if __name__ == "__main__":
    main()