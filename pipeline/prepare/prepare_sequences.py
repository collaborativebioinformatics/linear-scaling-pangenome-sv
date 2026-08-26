"""
prepare_sequences.py - Build chr21 smoke-interval multi-FASTA for PGGB.

Extracts the configured 0-based half-open interval from the reference
(via direct slicing — single contig) and from each HPRC haplotype
(via samtools faidx on the exact source_contig, using source_start/source_end
from the mapping report).

NEVER concatenates all contigs of a multi-contig assembly.
NEVER loads entire genomes into Python memory.
"""
import csv
import os
import sys
import yaml

from pipeline.prepare.faidx_utils import (
    faidx_extract, ensure_faidx, _revcomp, extract_whole_contig)


def find_assembly_path(name):
    """Find an assembly file by name in expected directories."""
    for d in ["work/downloads", "/data/hprc"]:
        for ext in [".fa.gz", ".fa", ".fasta", ".fna"]:
            p = os.path.join(d, name + ext)
            if os.path.exists(p):
                return p
    return None


def main():
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    tgt = cfg["target"]
    chrom, rs, re = tgt["chromosome"], tgt["start"], tgt["end"]
    print(f"=== Prepare Sequences ===")
    print(f"  Interval: {chrom}:{rs}-{re} ({re-rs} bp, 0-based half-open)")
    os.makedirs("results/preparation", exist_ok=True)
    os.makedirs("work/preparation", exist_ok=True)

    # Load mapping — must have exactly 4 mapped with source coordinates
    mp = "results/preparation/sequence_mapping.tsv"
    if not os.path.exists(mp):
        print("No mapping. Run: python3 pipeline/prepare/map_chromosome.py")
        sys.exit(1)
    mapping = [r for r in csv.DictReader(open(mp), delimiter="\t")
               if r["status"] == "mapped"]
    if len(mapping) != 4:
        print(f"FATAL: Expected 4 mapped, got {len(mapping)}.", file=sys.stderr)
        sys.exit(1)
    for r in mapping:
        if not r.get("source_start") or not r.get("source_end"):
            print(f"FATAL: {r['sample']} ({r['haplotype_label']}) "
                  f"missing source coordinates.", file=sys.stderr)
            sys.exit(1)

    # Reference — single contig, slice directly
    ref = "work/reference/GRCh38_chr21.fa"
    if not os.path.exists(ref):
        print("Ref not found."); sys.exit(1)

    multi = "results/preparation/chr21_multi.fa"
    pc = 0
    validation = []
    with open(multi, "w") as out:
        # Reference (path 1) — read the single contig and slice
        with open(ref) as f:
            ref_seq = "".join(line.strip() for line in f if not line.startswith(">"))
        ref_chunk = ref_seq[rs:re]
        out.write(f">GRCh38#0#{chrom}\n")
        for i in range(0, len(ref_chunk), 80):
            out.write(ref_chunk[i:i + 80] + "\n")
        pc += 1
        validation.append(("GRCh38#0#chr21", "GRCh38", rs, re, "+", len(ref_chunk)))
        print(f"  [1/5] GRCh38#0#{chrom} ({len(ref_chunk)} bp)")

        # HPRC haplotypes — extract from exact source_contig via samtools faidx
        for r in mapping:
            sm = r["sample"]
            nh = r["haplotype"]
            hl = r["haplotype_label"]
            name = r["assembly_name"]
            contig = r["source_contig"]
            ss = int(r["source_start"])
            se = int(r["source_end"])
            strand = r["strand"]

            ap = find_assembly_path(name)
            if not ap:
                print(f"FATAL: {name} not found", file=sys.stderr)
                sys.exit(1)

            # Ensure FASTA is indexed
            ensure_faidx(ap)

            # Extract from the exact source_contig using samtools faidx
            # PAF coordinates are 0-based half-open
            chunk = faidx_extract(ap, contig, ss, se, strand)

            out.write(f">{sm}#{nh}#{chrom}\n")
            for i in range(0, len(chunk), 80):
                out.write(chunk[i:i + 80] + "\n")
            pc += 1
            validation.append((f"{sm}#{nh}#{chrom}", contig, ss, se, strand, len(chunk)))
            print(f"  [{pc}/5] {sm}#{nh}#{chrom} ({hl}, "
                  f"{len(chunk)} bp from {contig}, strand={strand})")

    # Validation table
    print(f"\n{'Path':<30} {'Contig':<30} {'Strand':<8} {'Length':<10}")
    print(f"{'-'*30} {'-'*30} {'-'*8} {'-'*10}")
    for v in validation:
        print(f"{v[0]:<30} {v[1]:<30} {v[4]:<8} {v[5]:<10}")

    # Reject if any HPRC sequence >2x or <0.1x reference interval length
    ref_len = validation[0][5]
    for v in validation[1:]:
        if v[5] > ref_len * 2:
            print(f"\nFATAL: {v[0]} length ({v[5]} bp) >2x ref "
                  f"({ref_len} bp). Bad mapping.", file=sys.stderr)
            sys.exit(1)
        if v[5] < ref_len * 0.1:
            print(f"\nFATAL: {v[0]} length ({v[5]} bp) <0.1x ref "
                  f"({ref_len} bp). Bad mapping or partial coverage.",
                  file=sys.stderr)
            sys.exit(1)

    sz = os.path.getsize(multi) / 1e6
    print(f"\nOutput: {pc} paths, {sz:.1f} MB -> {multi}")
    if pc != 5:
        print("FATAL: Expected 5 paths.", file=sys.stderr)
        sys.exit(1)
    print("Exactly 5 paths. Ready for PGGB.")


if __name__ == "__main__":
    main()