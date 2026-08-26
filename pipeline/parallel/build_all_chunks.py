#!/usr/bin/env python3
"""
build_all_chunks.py — Build all chunk graphs, locally or in parallel.

TODO (Quang): Add subprocess worker pool for local parallel builds.
TODO (Ali): Add DNAnexus job submission mode.

Usage:
    python3 pipeline/parallel/build_all_chunks.py
"""

import csv
import os


def main():
    manifest_path = "work/chunks/chunk_manifest.tsv"
    if not os.path.exists(manifest_path):
        print(f"Chunk manifest not found: {manifest_path}")
        print("Run: python3 pipeline/parallel/make_chunks.py")
        return

    with open(manifest_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        chunks = list(reader)

    if not chunks:
        print("No chunks in manifest.")
        return

    print(f"Found {len(chunks)} chunks")

    for chunk in chunks:
        cid = chunk["chunk_id"]
        gfa_path = f"work/chunks/{cid}.gfa"
        if os.path.exists(gfa_path):
            print(f"  SKIP {cid}: already exists")
            continue
        print(f"  SUBMIT {cid}: chr21:{chunk['reference_start']}-{chunk['reference_end']}")
        print(f"    REAL-DATA STEP: bash pipeline/parallel/build_chunk.sh {cid}")


if __name__ == "__main__":
    main()