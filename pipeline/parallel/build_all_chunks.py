"""
build_all_chunks.py - Build per-chunk FASTA files with orthologous intervals.

For each chunk, slices the REFERENCE by coordinate, then maps the same
interval to each HPRC haplotype using the pre-computed mapping.
Every chunk FASTA contains exactly 5 LOCALLY corresponding sequences.
"""
import argparse
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


def load_manifest(path="work/chunks/chunk_manifest.tsv"):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_mapping(path="results/preparation/sequence_mapping.tsv"):
    if not os.path.exists(path):
        return []
    return list(csv.DictReader(open(path), delimiter="\t"))


def get_seq(ap):
    """Read entire assembly FASTA into a single string."""
    seq = []
    with _open_read(ap) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def map_chunk_to_source(rs, re, ss, se, cs, ce, strand):
    """Map a ref chunk interval [cs,ce) to source coordinates via linear scaling.

    All coordinates are 0-based half-open.
    Returns (src_start, src_end).
    """
    ref_len = re - rs
    if ref_len == 0:
        return (ss, ss)
    frac_s = (cs - rs) / ref_len
    frac_e = (ce - rs) / ref_len
    src_s = int(ss + frac_s * (se - ss))
    src_e = int(ss + frac_e * (se - ss))
    if strand == "-":
        orig_s, orig_e = src_s, src_e
        src_s = se - (orig_e - ss)
        src_e = se - (orig_s - ss)
    return (src_s, src_e)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true",
                        help="Run PGGB on each chunk")
    args = parser.parse_args()

    chunks = load_manifest()
    if not chunks:
        print("No chunk manifest. Run make_chunks.py first.")
        return

    mapping = load_mapping()
    if not mapping:
        print("No mapping. Run map_chromosome.py first.")
        return

    ref = "work/reference/GRCh38_chr21.fa"
    if not os.path.exists(ref):
        print("Ref not found.")
        return
    ref_seq = get_seq(ref)

    # Pre-load HPRC assembly sequences
    hap_seqs = {}
    for r in mapping:
        name = r["assembly_name"]
        if name in hap_seqs:
            continue
        ap = None
        for d in ["work/downloads", "/data/hprc"]:
            for ext in [".fa.gz", ".fa", ".fasta", ".fna"]:
                p = os.path.join(d, name + ext)
                if os.path.exists(p):
                    ap = p
                    break
            if ap:
                break
        if ap:
            hap_seqs[name] = get_seq(ap)
        else:
            print(f"FATAL: {name} not found")
            sys.exit(1)

    os.makedirs("work/chunks", exist_ok=True)
    print(f"Building {len(chunks)} orthologous chunk FASTA files...")

    for c in chunks:
        cid = c["chunk_id"]
        cs = int(c["reference_start"])
        ce = int(c["reference_end"])
        op = f"work/chunks/{cid}.fa"

        if os.path.exists(op) and os.path.getsize(op) > 0:
            print(f"  EXISTS {cid}")
            continue

        with open(op, "w") as fout:
            # Reference: slice by coordinate
            ref_chunk = ref_seq[cs:ce]
            fout.write(f">GRCh38#0#chr21\n")
            for i in range(0, len(ref_chunk), 80):
                fout.write(ref_chunk[i:i+80] + "\n")

            # HPRC: map chunk interval to source coordinates
            for r in mapping:
                sm = r["sample"]
                nh = r["haplotype"]
                name = r["assembly_name"]
                ss = int(r["source_start"])
                se = int(r["source_end"])
                rs = int(r["reference_start"])
                re = int(r["reference_end"])
                strand = r["strand"]

                src_s, src_e = map_chunk_to_source(
                    rs, re, ss, se, cs, ce, strand)

                seq = hap_seqs[name][src_s:src_e]
                if strand == "-":
                    seq = _revcomp(seq)

                fout.write(f">{sm}#{nh}#chr21\n")
                for i in range(0, len(seq), 80):
                    fout.write(seq[i:i+80] + "\n")

        print(f"  {cid}: {cs}-{ce} -> {op}")

    if args.execute:
        print("\nBuilding chunk graphs via Docker...")
        for c in chunks:
            cid = c["chunk_id"]
            if not os.path.exists(f"work/chunks/{cid}.gfa"):
                subprocess.run(
                    ["bash", "pipeline/parallel/build_chunk.sh", cid],
                    check=False)

    print(f"Done. {len(chunks)} chunks.")


if __name__ == "__main__":
    main()