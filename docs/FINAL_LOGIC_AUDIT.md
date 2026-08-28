# FINAL LOGIC AUDIT — Static Results & Real 1 Mb Validation

## A. BASELINE 1 Mb VALIDATION (P0)

**Status: INVALID_FOR_1MB_BENCHMARK**

The current `results/baseline/baseline.gfa` was built from **full
chromosome-21 contigs**, not the intended 1 Mb slice.

| Sample#Hap | Spelled length | Steps | Subrange? |
|------------|---------------|-------|-----------|
| GRCh38#0 | 46,709,983 bp | 1 | No |
| HG00673#1 | 40,759,410 bp | 15,906 | No |
| HG00673#2 | 36,873,966 bp | 1 | No |
| HG00733#1 | 39,274,658 bp | 8,262 | No |
| HG00733#2 | 39,318,612 bp | 10,765 | No |

### WHY HG00673 PREVIOUSLY SHOWED ~40.7 MB

The path `HG00673#1#JAHBBZ020000061.1` is the **full HPRC contig**
(~40.7 Mb). The path name carries **no `:start-end` subrange**, which is
what PGGB emits when the input FASTA was sliced to a region. Therefore the
whole contig was fed to PGGB.

This is **not** a parser bug and **not** double-counting. It is genuinely
the full-contig sequence length. GRCh38 is a single 46.7 Mb node (the full
chromosome 21, and it is all `N`s — the reference was not even properly
extracted).

### Consequence

The 1,971-second timing in `run_metadata.json` is for a **full chr21**
build (~200 Mb across 5 paths), NOT a 1 Mb build. It must NOT be used as
the official 1 Mb baseline. The correct baseline must be rebuilt.

## B. FINAL BENCHMARK REGION

- **GRCh38 chr21:20,000,000-21,000,000** (1,000,000 bp)
- 5 haplotypes: GRCh38, HG00673 (hap1/2), HG00733 (hap1/2)

## C. FINAL VISUALIZATION REGION

Not yet selected for the real 1 Mb benchmark (the 1 Mb baseline does not
exist yet). `pipeline/export/select_visualization_region.py` is ready to
pick a 50 kb high-diversity window once the real 1 Mb graph exists.

## D. TIMING

| Metric | Value |
|--------|-------|
| Baseline wall time | NOT_RUN (current 1,971s is full-chr21, invalid) |
| Parallel wall time | NOT_RUN |
| Speedup | NOT_RUN |
| Sum worker time | NOT_RUN |
| Orchestration | NOT_RUN |

## E. GRAPH / PATH / BOUNDARY COMPARISON

| Comparison | Status |
|-----------|--------|
| Baseline nodes/edges | full-chr21 only (invalid) |
| Stitched nodes/edges | placeholder (synthetic chunks, linear only) |
| Exact path matches | NOT_RUN |
| Boundary status | NOT_RUN |

## F. TRUVARI

| Comparison | TP | FP | FN | Prec | Recall | F1 | Status |
|-----------|----|----|----|------|--------|----|--------|
| baseline vs stitched | — | — | — | — | — | — | NOT_RUN |
| baseline vs dipcall | — | — | — | — | — | — | NOT_RUN |
| stitched vs dipcall | — | — | — | — | — | — | NOT_RUN |

## G. VALIDATION LEVEL

| Item | Level |
|------|-------|
| Algorithm | IMPLEMENTED |
| Synthetic | SYNTHETIC_VALIDATED (demo) |
| Real HPRC baseline | INVALID_FOR_1MB_BENCHMARK (full contigs) |
| Real HPRC parallel | NOT_RUN |
| Real stitch | NOT_RUN |
| Real graph equivalence | NOT_RUN |
| Real variant equivalence | NOT_RUN |

## H. WEB BEHAVIOR

- Static data only: YES
- DNAnexus required at view time: NO
- Recompute on sample change: NO
- Recompute on build: NO

## I. DNANEXUS

- Auth: EXPIRED (token timeout). Requires `dx login`.
- Required project: `project-JB6zQBj0ZQv2Bk79ggBBv76Z`
- Next steps (gated): build applets → dry-run → micro-test → real 1 Mb run.

## J. WHAT STILL REQUIRES DNANEXUS

1. Rebuild the **real 1 Mb baseline** (currently full contigs).
2. Run **parallel chunks** with canonical config (`-K 19`, not `-k 19`).
3. Run **stitch** on real chunk GFAs.
4. Run **vg deconstruct + Truvari** on the real 1 Mb graphs.
5. Record **timing** from real DNAnexus job metadata.
