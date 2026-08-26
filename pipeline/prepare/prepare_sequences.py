"""
prepare_sequences.py - Build chr21 smoke-interval multi-FASTA for PGGB.

Extracts the configured 0-based half-open interval from both reference
and each HPRC haplotype (using source_start/source_end from mapping).
Reverse-complements sequences mapped on the negative strand.
Requires exactly 5 paths. Fails otherwise.
"""
import csv
import gzip
import os
import sys
import yaml


def _revcomp(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G",
            "a": "t", "t": "a", "g": "c", "c": "g",
            "N": "N", "n": "n"}
    return "".join(comp.get(c, c) for c in reversed(seq))


def _open_read(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def extract_region(fa, start, end, strand="+"):
    """Extract a 0-based half-open interval, RC if strand=-."""
    seq = []
    with _open_read(fa) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    s = "".join(seq)
    chunk = s[start:end]
    if strand == "-":
        chunk = _revcomp(chunk)
    return chunk


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

    ref = "work/reference/GRCh38_chr21.fa"
    if not os.path.exists(ref):
        print("Ref not found."); sys.exit(1)

    multi = "results/preparation/chr21_multi.fa"
    pc = 0
    validation = []
    with open(multi, "w") as out:
        # Reference (path 1) — slice the interval
        ref_chunk = extract_region(ref, rs, re)
        out.write(f">GRCh38#0#{chrom}\n")
        for i in range(0, len(ref_chunk), 80):
            out.write(ref_chunk[i:i+80] + "\n")
        pc += 1
        validation.append(("GRCh38#0#chr21", rs, re, "GRCh38",
                          rs, re, "+", len(ref_chunk)))
        print(f"  [1/5] GRCh38#0#{chrom} ({len(ref_chunk)} bp)")

        for r in mapping:
            sm = r["sample"]
            nh = r["haplotype"]
            hl = r["haplotype_label"]
            name = r["assembly_name"]
            ss = int(r["source_start"])
            se = int(r["source_end"])
            strand = r["strand"]

            # Find assembly file
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
                print(f"FATAL: {name} not found", file=sys.stderr)
                sys.exit(1)

            # Extract interval, RC if strand==-
            chunk = extract_region(ap, ss, se, strand)
            out.write(f">{sm}#{nh}#{chrom}\n")
            for i in range(0, len(chunk), 80):
                out.write(chunk[i:i+80] + "\n")
            pc += 1
            validation.append((f"{sm}#{nh}#{chrom}", ss, se, name,
                              ss, se, strand, len(chunk)))
            print(f"  [{pc}/5] {sm}#{nh}#{chrom} ({hl}, "
                  f"{len(chunk)} bp, strand={strand})")

    # Validation table
    print(f"\n{'Path':<30} {'Interval':<28} {'Strand':<8} {'Length':<10}")
    print(f"{'-'*30} {'-'*28} {'-'*8} {'-'*10}")
    for v in validation:
        tint = f"{chrom}:{v[1]}-{v[2]}"
        print(f"{v[0]:<30} {tint:<28} {v[6]:<8} {v[7]:<10}")

    # Reject if any HPRC sequence >2x reference interval
    ref_len = validation[0][7]
    for v in validation[1:]:
        if v[7] > ref_len * 2:
            print(f"\nFATAL: {v[0]} length ({v[7]} bp) >2x ref "
                  f"({ref_len} bp). Bad mapping.", file=sys.stderr)
            sys.exit(1)

    sz = os.path.getsize(multi) / 1e6
    print(f"\nOutput: {pc} paths, {sz:.1f} MB -> {multi}")
    if pc != 5:
        print("FATAL: Expected 5 paths.", file=sys.stderr)
        sys.exit(1)
    print("Exactly 5 paths. Ready for PGGB.")


if __name__ == "__main__":
    main()