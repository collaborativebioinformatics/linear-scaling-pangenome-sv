# Zero-Edge Chunk Root Cause

## SYMPTOM

All 3 real 1 Mb chunks produced GFAs with 5 segments, 5 paths, and **0 edges**
regardless of using `mash_kmer=31` or `mash_kmer=19`. The monolithic baseline
(10,389 nodes, 14,185 edges) built correctly from the same container.

## ROOT CAUSE

PGGB's `pggb` wrapper uses a **two-pass wfmash** strategy:

1. **Pass 1 (`--approx-map`)**: rapid mash-based mapping. Produces PAF records
   with **no CIGAR** — coordinates are approximate, no exact alignment.

2. **Pass 2 (`--lower-triangular --invert-filtering`)**: full base-level
   alignment only on haplotype pairs **NOT already covered** by pass 1.
   Produces PAF records with **detailed CIGAR**.

The two PAFs are concatenated and fed to seqwish. Seqwish needs CIGAR strings
to find exact matches (controlled by `-k 29`). Without CIGAR, seqwish cannot
find any exact 29-mer matches across segment boundaries.

### Why Baseline Worked (5 Mb, 5 haplotypes)

With 5,006,733 bp across 5 haplotypes, pass 1's mash mapping does NOT cover
every pair exhaustively. Pass 2 fills in the remaining pairs with detailed
alignments. The combined PAF has both mapping records AND detailed alignment
records with CIGAR. Seqwish finds 12,448 links.

### Why Chunks Failed (2.1 Mb, 5 haplotypes)

With only 2,127,860 bp across 5 haplotypes, pass 1's mash mapping covers
ALL pairs. Pass 2 finds NOTHING to do (every pair was already mapped in
pass 1). The combined PAF contains only mash-based records with NO CIGAR.
Seqwish finds 0 links.

### Confirmed in Logs

```
Baseline: [seqwish::links] links derived (12448 links)
Chunk:    [seqwish::links] links derived (0 links)
```

Both used identical wfmash commands — the difference is purely the ratio
of input size to mash coverage, making pass 2 unnecessary (and therefore
absent) for small chunks.

## FIX

Use a **single-pass detailed wfmash** (with `-L` for lower-triangular, no
`--approx-map`, no `--invert-filtering`), producing full base-level
alignments with CIGAR in one pass. Feed this directly to seqwish.

Smoothxg flags must exactly match the validated baseline:
```
smoothxg -t 8 -T 8 -g seqwish.gfa -r 5 --base ... --chop-to 100
  -I .9000 -R 0 -j 0 -e 0 -l 700,1100 -p 1,4,6,2,26,1
  -O 0.001 -Y 500 -d 0 -D 0 -Q Consensus_ -V -o smooth.gfa
```

The previous manual pipeline had two bugs:
- `smoothxg -m /data/input.fa` — `-m` means "write MAF", not "input FASTA"
- `smoothxg -n 5` — should be `-r 5` (haplotype count)

## CHUNK APPLET ID

`applet-JB8qQj80ZQvPyVkGzqxb5bVY`