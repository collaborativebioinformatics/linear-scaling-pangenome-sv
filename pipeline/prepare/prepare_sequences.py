#!/usr/bin/env python3
"""
prepare_sequences.py - Build the smoke-interval multi-FASTA for PGGB.

Extracts the configured 0-based half-open interval from the reference and
from each HPRC haplotype, using samtools faidx against the exact
source_contig / source_start / source_end recorded by map_chromosome.py.

Guarantees:
  - NEVER concatenates all contigs of a multi-contig assembly.
  - NEVER loads an entire genome into Python memory (the reference is
    extracted with faidx too, not slurped and sliced).
  - Emits PanSN-named paths: sample#haplotype#contig.

Output name encodes the interval so different intervals cannot silently
overwrite one another, and a stable `<chrom>_multi.fa` symlink/copy is
also written for tools that expect a fixed name.

Usage:
    python3 pipeline/prepare/prepare_sequences.py
    python3 pipeline/prepare/prepare_sequences.py --config config/pipeline.yaml
    python3 pipeline/prepare/prepare_sequences.py --reference work/reference/GRCh38_chr21.fa
"""
import argparse
import csv
import os
import shutil
import sys

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from pipeline.prepare.faidx_utils import (  # noqa: E402
    ensure_faidx, faidx_extract, get_contig_length)

LINE_WIDTH = 80
SEARCH_DIRS = ["work/downloads", "/data/hprc"]
FASTA_EXTS = [".fa.gz", ".fa", ".fasta", ".fna"]


def find_assembly_path(name):
    """Find an assembly file by assembly_name in the expected directories."""
    for d in SEARCH_DIRS:
        for ext in FASTA_EXTS:
            p = os.path.join(d, name + ext)
            if os.path.exists(p):
                return p
    return None


def reference_name(cfg):
    """Reference assembly name from config, e.g. "GRCh38".

    This is `target.reference` in config/pipeline.yaml — the same value
    scripts/prepare_reference.sh uses to name its output and the same
    string that becomes the PanSN sample field of the reference path.
    """
    return (cfg.get("target") or {}).get("reference", "GRCh38")


def find_reference(cfg, chrom, explicit=None):
    """Locate the reference FASTA for this chromosome.

    Order: --reference flag, then the project convention
    work/reference/<target.reference>_<chrom>.fa (what
    scripts/prepare_reference.sh writes), then any single FASTA in
    work/reference/. Both the reference name and the chromosome come from
    config, so switching target.reference to CHM13 or target.chromosome to
    chr6 keeps working.
    """
    if explicit:
        return explicit
    conventional = "work/reference/%s_%s.fa" % (reference_name(cfg), chrom)
    if os.path.exists(conventional):
        return conventional
    d = "work/reference"
    if os.path.isdir(d):
        fas = sorted(f for f in os.listdir(d)
                     if f.endswith((".fa", ".fasta", ".fa.gz")))
        if len(fas) == 1:
            return os.path.join(d, fas[0])
    return conventional  # report this path in the error message


def write_fasta_record(handle, name, seq, width=LINE_WIDTH):
    handle.write(">%s\n" % name)
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def resolve_reference_contig(ref_path, chrom):
    """Pick the contig in the reference matching the target chromosome.

    Accepts `chr21`, `21`, or a PanSN-prefixed `GRCh38#0#chr21`.
    """
    fai = ref_path + ".fai"
    ensure_faidx(ref_path)
    names = []
    with open(fai) as f:
        for line in f:
            if line.strip():
                names.append(line.split("\t")[0])
    if not names:
        raise RuntimeError("empty .fai for %s" % ref_path)
    for cand in (chrom, chrom.replace("chr", ""), "chr" + chrom.lstrip("chr")):
        if cand in names:
            return cand
    for n in names:                       # PanSN-prefixed reference
        if n.endswith("#" + chrom) or n.endswith(chrom):
            return n
    if len(names) == 1:
        return names[0]                   # single-contig reference file
    raise RuntimeError(
        "could not find contig %r in %s (contigs: %s)"
        % (chrom, ref_path, ", ".join(names[:5])))


def main():
    ap = argparse.ArgumentParser(
        description="Build the interval multi-FASTA for PGGB.")
    ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--reference", default=None,
                    help="Reference FASTA (overrides config and convention).")
    ap.add_argument("--mapping",
                    default="results/preparation/sequence_mapping.tsv")
    ap.add_argument("--outdir", default="results/preparation")
    ap.add_argument("--max-ratio", type=float, default=2.0,
                    help="Reject a haplotype longer than this x reference.")
    ap.add_argument("--min-ratio", type=float, default=0.1,
                    help="Reject a haplotype shorter than this x reference.")
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    tgt = cfg["target"]
    chrom = tgt["chromosome"]
    rs, re_ = int(tgt["start"]), int(tgt["end"])

    print("=== Prepare Sequences ===")
    print("  Interval: %s:%d-%d (%d bp, 0-based half-open)"
          % (chrom, rs, re_, re_ - rs))

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs("work/preparation", exist_ok=True)

    # --- mapping report -----------------------------------------------------
    if not os.path.exists(args.mapping):
        print("No mapping. Run: python3 pipeline/prepare/map_chromosome.py",
              file=sys.stderr)
        sys.exit(1)
    with open(args.mapping) as f:
        mapping = [r for r in csv.DictReader(f, delimiter="\t")
                   if r.get("status") == "mapped"]
    if not mapping:
        print("FATAL: no mapped haplotypes in %s" % args.mapping,
              file=sys.stderr)
        sys.exit(1)
    for r in mapping:
        if not r.get("source_start") or not r.get("source_end"):
            print("FATAL: %s (%s) missing source coordinates."
                  % (r.get("sample"), r.get("haplotype_label")),
                  file=sys.stderr)
            sys.exit(1)

    expected_paths = len(mapping) + 1     # haplotypes + reference
    print("  Haplotypes mapped: %d (expecting %d paths incl. reference)"
          % (len(mapping), expected_paths))

    # --- reference ----------------------------------------------------------
    ref = find_reference(cfg, chrom, args.reference)
    if not os.path.exists(ref):
        print("FATAL: reference not found: %s" % ref, file=sys.stderr)
        print("  Run: bash scripts/prepare_reference.sh", file=sys.stderr)
        sys.exit(1)
    ensure_faidx(ref)
    ref_contig = resolve_reference_contig(ref, chrom)
    ref_len_total = get_contig_length(ref, ref_contig)
    if re_ > ref_len_total:
        print("FATAL: interval end %d exceeds %s length %d"
              % (re_, ref_contig, ref_len_total), file=sys.stderr)
        sys.exit(1)

    # Interval-stamped name so two intervals cannot overwrite each other.
    stamped = os.path.join(args.outdir,
                           "%s_%d_%d_multi.fa" % (chrom, rs, re_))
    stable = os.path.join(args.outdir, "%s_multi.fa" % chrom)

    validation = []
    with open(stamped, "w") as out:
        # Reference via faidx — do NOT slurp the file into memory.
        ref_chunk = faidx_extract(ref, ref_contig, rs, re_, "+")
        ref_name = "%s#0#%s" % (reference_name(cfg), chrom)
        write_fasta_record(out, ref_name, ref_chunk)
        validation.append((ref_name, ref_contig, rs, re_, "+", len(ref_chunk)))
        print("  [1/%d] %s (%d bp)" % (expected_paths, ref_name,
                                       len(ref_chunk)))

        for i, r in enumerate(mapping, start=2):
            sm = r["sample"]
            nh = r["haplotype"]
            hl = r.get("haplotype_label", "")
            name = r["assembly_name"]
            contig = r["source_contig"]
            ss, se = int(r["source_start"]), int(r["source_end"])
            strand = r.get("strand", "+")

            ap_path = find_assembly_path(name)
            if not ap_path:
                print("FATAL: assembly not found: %s" % name, file=sys.stderr)
                print("  Searched: %s" % ", ".join(SEARCH_DIRS),
                      file=sys.stderr)
                print("  Run: python3 scripts/download_hprc.py --execute",
                      file=sys.stderr)
                sys.exit(1)

            ensure_faidx(ap_path)
            chunk = faidx_extract(ap_path, contig, ss, se, strand)

            path_name = "%s#%s#%s" % (sm, nh, chrom)
            write_fasta_record(out, path_name, chunk)
            validation.append((path_name, contig, ss, se, strand, len(chunk)))
            print("  [%d/%d] %s (%s, %d bp from %s, strand=%s)"
                  % (i, expected_paths, path_name, hl, len(chunk), contig,
                     strand))

    # --- validation ---------------------------------------------------------
    print("\n%-30s %-30s %-8s %-10s"
          % ("Path", "Contig", "Strand", "Length"))
    print("%s %s %s %s" % ("-" * 30, "-" * 30, "-" * 8, "-" * 10))
    for v in validation:
        print("%-30s %-30s %-8s %-10s" % (v[0], v[1], v[4], v[5]))

    ref_len = validation[0][5]
    if ref_len == 0:
        print("\nFATAL: reference interval extracted 0 bp.", file=sys.stderr)
        sys.exit(1)
    for v in validation[1:]:
        ratio = v[5] / ref_len
        if ratio > args.max_ratio:
            print("\nFATAL: %s is %.2fx the reference interval (%d vs %d bp). "
                  "Bad mapping." % (v[0], ratio, v[5], ref_len),
                  file=sys.stderr)
            sys.exit(1)
        if ratio < args.min_ratio:
            print("\nFATAL: %s is %.2fx the reference interval (%d vs %d bp). "
                  "Bad mapping or partial coverage."
                  % (v[0], ratio, v[5], ref_len), file=sys.stderr)
            sys.exit(1)

    n_paths = len(validation)
    if n_paths != expected_paths:
        print("FATAL: wrote %d paths, expected %d."
              % (n_paths, expected_paths), file=sys.stderr)
        sys.exit(1)

    # Stable alias so downstream steps that hardcode <chrom>_multi.fa work.
    # build_baseline.sh defaults differ between branches — one expects
    # chr21_multi.fa, another chr21_20000000_21000000_multi.fa. Writing both
    # means `make baseline` succeeds either way.
    if os.path.abspath(stable) != os.path.abspath(stamped):
        shutil.copyfile(stamped, stable)

    sz = os.path.getsize(stamped) / 1e6
    print("\nOutput: %d paths, %.1f MB" % (n_paths, sz))
    print("  %s" % stamped)
    print("  %s  (stable alias)" % stable)
    print("Ready for PGGB.")


if __name__ == "__main__":
    main()
