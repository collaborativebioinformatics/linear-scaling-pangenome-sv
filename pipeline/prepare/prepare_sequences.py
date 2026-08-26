"""
prepare_sequences.py — Build chr21 multi-FASTA for PGGB.

NO naive coordinate slicing. Uses map_chromosome.py to identify correct
chr21 contigs, extracts them fully. PGGB handles alignment.
"""
import csv
import os
import sys
import yaml


def main():
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    chrom = cfg["target"]["chromosome"]
    print(f"=== Prepare Sequences (target: {chrom}) ===")
    os.makedirs("results/preparation", exist_ok=True)
    os.makedirs("work/preparation", exist_ok=True)

    mp = "results/preparation/sequence_mapping.tsv"
    if not os.path.exists(mp):
        print("No mapping. Run: python3 pipeline/prepare/map_chromosome.py")
        sys.exit(1)
    mapping = list(csv.DictReader(open(mp), delimiter="\t"))

    ref = "work/reference/GRCh38_chr21.fa"
    if not os.path.exists(ref):
        print("Reference not found. Run: bash scripts/prepare_reference.sh")
        sys.exit(1)

    multi = "results/preparation/chr21_multi.fa"
    pc = 0
    with open(multi, "w") as out:
        with open(ref) as f:
            out.write(f.read())
        pc += 1
        print(f"  GRCh38#0#chr21 (ref)")

        for row in mapping:
            if row["status"] != "mapped":
                print(f"  SKIP {row['sample']} ({row['haplotype']}) — not mapped")
                continue

            sm, hap = row["sample"], row["haplotype"]
            tag = "mat" if hap == "maternal" else "pat"
            name_key = f"{sm}_{tag}_hprc_r2_v1.0.1"
            contig_name = row["source_contig"].split()[0]

            ap = None
            for d in ["work/downloads", "/data/hprc"]:
                for ext in ["", ".fa", ".fasta", ".fna", ".fa.gz"]:
                    p = os.path.join(d, name_key + ext)
                    if os.path.exists(p):
                        ap = p
                        break
                if ap:
                    break
            if not ap:
                print(f"  SKIP {sm} ({hap}) — file not found")
                continue

            seq_path = f"work/preparation/{sm}_{hap}_chr21.fa"
            with open(ap) as fin, open(seq_path, "w") as fout:
                mode = False
                for line in fin:
                    if line.startswith(">"):
                        h = line[1:].strip().split()[0]
                        if (h == contig_name or h.endswith("chr21")
                                or "#chr21" in h):
                            mode = True
                            fout.write(f">{sm}#{hap}#{chrom}\n")
                        else:
                            mode = False
                    elif mode:
                        fout.write(line)

            if os.path.getsize(seq_path) > 50:
                with open(seq_path) as sf:
                    next(sf)  # skip header
                    for line in sf:
                        out.write(line)
                pc += 1
                sz = os.path.getsize(seq_path) / 1e6
                print(f"  {sm}#{hap}#{chrom} ({sz:.1f} MB)")
            else:
                # Fallback: write full assembly
                print(f"  WARNING: chr21 not found in {sm} ({hap}). "
                      "Writing full assembly.")
                with open(ap) as fin:
                    out.write(f">{sm}#{hap}#{chrom}\n")
                    for line in fin:
                        if not line.startswith(">"):
                            out.write(line)
                pc += 1

    sz = os.path.getsize(multi) / 1e6
    print(f"\nOutput: {pc} paths, {sz:.1f} MB -> {multi}")
    print("NOTE: No naive coordinate slicing. PGGB handles alignment.")
    print(f"Next: bash pipeline/baseline/build_baseline.sh {multi}")


if __name__ == "__main__":
    main()