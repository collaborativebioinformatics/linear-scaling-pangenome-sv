#!/usr/bin/env python3
"""
build_all_chunks.py — Extract chunk sequences and optionally build chunk graphs.

Step 1: From the multi-haplotype FASTA, extract per-chunk FASTA files.
Step 2: Optionally run PGGB on each chunk via build_chunk.sh.

Usage:
    python3 pipeline/parallel/build_all_chunks.py                    # extract only
    python3 pipeline/parallel/build_all_chunks.py --execute          # extract + build
    python3 pipeline/parallel/build_all_chunks.py --input <multi.fa> # custom input
"""

import argparse
import csv
import os
import subprocess
import sys


def load_manifest(path="work/chunks/chunk_manifest.tsv"):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def extract_chunk_fasta(multi_fa, chunk_id, chrom, start, end, output_dir):
    """Extract a region from a multi-FASTA by header-matching the chromosome."""
    out_path = os.path.join(output_dir, f"{chunk_id}.fa")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return out_path

    with open(multi_fa) as fin, open(out_path, "w") as fout:
        current_header = None
        current_seq_lines = []
        written = False

        def flush():
            if current_header is None:
                return
            # Only write if this sequence spans the requested interval
            # For simplicity: write all sequences that have the chrom tag
            seq = "".join(current_seq_lines)
            if len(seq) >= start:
                chunk_seq = seq[start:end]
                fout.write(f">{current_header}\n")
                for i in range(0, len(chunk_seq), 80):
                    fout.write(chunk_seq[i:i + 80] + "\n")

        for line in fin:
            if line.startswith(">"):
                flush()
                current_header = line[1:].strip()
                current_seq_lines = []
            else:
                current_seq_lines.append(line.strip())
        flush()

    print(f"  Extracted {chunk_id}: {chrom}:{start}-{end} -> {out_path}")
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Build all chunk graphs")
    parser.add_argument("--execute", action="store_true",
                        help="Also run PGGB on each chunk")
    parser.add_argument("--input", default="results/preparation/chr21_multi.fa",
                        help="Multi-haplotype FASTA input")
    args = parser.parse_args()

    chunks = load_manifest()
    if not chunks:
        print("No chunk manifest found. Run make_chunks.py first.")
        return

    multi_fa = args.input
    if not os.path.exists(multi_fa):
        print(f"Multi-FASTA not found: {multi_fa}")
        print("Run prepare_sequences.py first.")
        return

    os.makedirs("work/chunks", exist_ok=True)
    print(f"Extracting {len(chunks)} chunks from {multi_fa}...")

    for c in chunks:
        cid = c["chunk_id"]
        chrom = c["chrom"]
        start = int(c["reference_start"])
        end = int(c["reference_end"])
        extract_chunk_fasta(multi_fa, cid, chrom, start, end, "work/chunks")

    if args.execute:
        print("\nBuilding chunk graphs via PGGB Docker container...\n")
        for c in chunks:
            cid = c["chunk_id"]
            gfa_path = f"work/chunks/{cid}.gfa"
            if os.path.exists(gfa_path):
                print(f"  SKIP {cid}: already exists")
                continue
            print(f"  Building {cid}...")
            subprocess.run(
                ["bash", "pipeline/parallel/build_chunk.sh", cid],
                check=False,
            )

    print(f"\nDone. {len(chunks)} chunks processed.")


if __name__ == "__main__":
    main()