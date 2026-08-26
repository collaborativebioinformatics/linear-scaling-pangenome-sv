"""
prepare_sequences.py — Build chr21 multi-FASTA for PGGB.

NO naive coordinate slicing. Uses map_chromosome.py to identify correct
chr21 contigs, extracts them fully. Uses gzip.open for compressed assemblies.

Requires exactly 5 paths (GRCh38 + 4 HPRC haplotypes). Fails otherwise.
"""
import csv
import gzip
import os
import sys
import yaml


def _open_read(path):
    """Open a FASTA file for reading, handling .gz transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main():
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    chrom = cfg["target"]["chromosome"]
    print(f"=== Prepare Sequences (target: {chrom}) ===")
    os.makedirs("results/preparation", exist_ok=True)
    os.makedirs("work/preparation", exist_ok=True)

    # Load mapping — must have exactly 4 mapped haplotypes
    mp = "results/preparation/sequence_mapping.tsv"
    if not os.path.exists(mp):
        print("No mapping. Run: python3 pipeline/prepare/map_chromosome.py")
        sys.exit(1)
    mapping = [r for r in csv.DictReader(open(mp), delimiter="\t")
               if r["status"] == "mapped"]

    if len(mapping) != 4:
        print(f"FATAL: Expected 4 mapped haplotypes, got {len(mapping)}.",
              file=sys.stderr)
        print("All 4 HPRC haplotypes must be mapped to chr21. "
              "Check sequence_mapping.tsv.", file=sys.stderr)
        sys.exit(1)

    # Reference must exist
    ref = "work/reference/GRCh38_chr21.fa"
    if not os.path.exists(ref):
        print("Reference not found. Run: bash scripts/prepare_reference.sh")
        sys.exit(1)

    multi = "results/preparation/chr21_multi.fa"
    path_count = 0

    with open(multi, "w") as out:
        # Reference (path 1)
        with open(ref) as f:
            out.write(f.read())
        path_count += 1
        print(f"  [1/5] GRCh38#0#chr21 (reference)")

        for row in mapping:
            sm, hap = row["sample"], row["haplotype"]
            tag = "mat" if hap == "maternal" else "pat"
            name_key = f"{sm}_{tag}_hprc_r2_v1.0.1"
            contig_name = row["source_contig"].split()[0]

            # Find assembly file (.fa.gz preferred, then .fa)
            ap = None
            for d in ["work/downloads", "/data/hprc"]:
                for ext in [".fa.gz", ".fa", ".fasta", ".fna"]:
                    p = os.path.join(d, name_key + ext)
                    if os.path.exists(p):
                        ap = p
                        break
                if ap:
                    break

            if not ap:
                print(f"  FATAL: Assembly file for {sm} ({hap}) not found "
                      f"at any expected path.", file=sys.stderr)
                sys.exit(1)

            # Extract the chr21 contig by matching header from mapping
            seq_path = f"work/preparation/{sm}_{hap}_chr21.fa"
            found = False
            with _open_read(ap) as fin, open(seq_path, "w") as fout:
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

            # Verify extraction produced something
            if os.path.getsize(seq_path) > 50:
                with open(seq_path) as sf:
                    next(sf)  # skip header (already written)
                    for line in sf:
                        out.write(line)
                path_count += 1
                sz = os.path.getsize(seq_path) / 1e6
                print(f"  [{path_count}/5] {sm}#{hap}#{chrom} ({sz:.1f} MB)")
            else:
                print(f"  FATAL: chr21 contig not found in {sm} ({hap}) "
                      f"assembly {name_key}. Cannot extract chr21 sequence.",
                      file=sys.stderr)
                print(f"  Checked contig name: {contig_name}", file=sys.stderr)
                print(f"  Assembly file: {ap}", file=sys.stderr)
                sys.exit(1)

    # Verify exactly 5 paths
    sz = os.path.getsize(multi) / 1e6
    print(f"\nOutput: {path_count} paths, {sz:.1f} MB -> {multi}")

    if path_count != 5:
        print(f"FATAL: Expected exactly 5 paths (GRCh38 + 4 haplotypes), "
              f"but got {path_count}.", file=sys.stderr)
        sys.exit(1)

    print("Exactly 5 paths verified. Ready for PGGB.")
    print(f"Next: bash pipeline/baseline/build_baseline.sh {multi}")


if __name__ == "__main__":
    main()