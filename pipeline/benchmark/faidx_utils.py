"""
faidx_utils.py — Extract intervals from indexed FASTA files using samtools.

PAF coordinates from minimap2 are 0-based half-open [start, end).
samtools faidx requires 1-based inclusive region format: contig:start-end.

Coordinate conversion:
    samtools_start = paf_start + 1
    samtools_end   = paf_end

This module NEVER loads entire assemblies into Python memory.
"""
import os
import subprocess
import sys
import tempfile


def _revcomp(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G",
            "a": "t", "t": "a", "g": "c", "c": "g",
            "N": "N", "n": "n"}
    return "".join(comp.get(c, c) for c in reversed(seq))


def ensure_faidx(fasta_path: str) -> None:
    """Ensure .fai and .gzi exist for a FASTA file."""
    fai = fasta_path + ".fai"
    gzi = fasta_path + ".gzi"
    if not os.path.exists(fai):
        print(f"    Building index: {fai}")
        subprocess.run(["samtools", "faidx", fasta_path],
                       check=True, capture_output=True)


def faidx_extract(fasta_path: str, contig: str,
                  start_0based: int, end_0based: int,
                  strand: str = "+") -> str:
    """Extract a 0-based half-open interval from an indexed FASTA via samtools.

    PAF coordinates are 0-based half-open [start, end).
    samtools faidx region = contig:start+1-end (1-based inclusive).

    Returns the extracted sequence, reverse-complemented if strand='-'.
    """
    samtools_start = start_0based + 1
    samtools_end = end_0based
    region = f"{contig}:{samtools_start}-{samtools_end}"

    try:
        result = subprocess.run(
            ["samtools", "faidx", fasta_path, region],
            capture_output=True, text=True, timeout=300)
    except FileNotFoundError:
        print("FATAL: samtools not found. Install with: "
              "sudo apt-get install -y samtools", file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print(f"FATAL: samtools faidx timed out on {region}", file=sys.stderr)
        sys.exit(1)

    if result.returncode != 0:
        print(f"FATAL: samtools faidx failed for {region} on {fasta_path}",
              file=sys.stderr)
        print(f"  stderr: {result.stderr.strip()}", file=sys.stderr)
        sys.exit(1)

    # Parse output: header line followed by sequence lines
    seq_parts = []
    for line in result.stdout.strip().split("\n"):
        if line.startswith(">"):
            continue
        seq_parts.append(line.strip())
    seq = "".join(seq_parts)

    if strand == "-":
        seq = _revcomp(seq)

    return seq


def extract_whole_contig(fasta_path: str, contig: str) -> str:
    """Extract an entire contig from an indexed FASTA."""
    return faidx_extract(fasta_path, contig, 0, 2**63 - 1)


def get_contig_length(fasta_path: str, contig: str) -> int:
    """Get contig length from .fai file without loading sequence."""
    fai_path = fasta_path + ".fai"
    if not os.path.exists(fai_path):
        ensure_faidx(fasta_path)
        fai_path = fasta_path + ".fai"

    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == contig:
                return int(parts[1])
    print(f"FATAL: Contig '{contig}' not found in {fai_path}",
          file=sys.stderr)
    sys.exit(1)