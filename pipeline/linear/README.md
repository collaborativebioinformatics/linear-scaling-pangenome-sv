# Linear / Assembly-Based SV Pipeline

Orthogonal validation for the pangenome graph pipeline. Calls structural
variants by aligning de novo HPRC assemblies directly against GRCh38,
with no graph in the loop.

## Why this exists

The core experiment compares two graphs:

```
Monolithic PGGB (baseline)  vs.  Parallel PGGB + Reassembly (merged)
```

That comparison is *internally* consistent — but if both graphs share a
construction artifact, comparing them to each other will never reveal it.
The assembly-based workflow supplies an independent opinion: variants
derived from minimap2 assembly alignments, using none of the PGGB
machinery. A variant found by the graph AND by dipcall is almost
certainly real. A variant found only by the graph deserves scrutiny.

**This is not the baseline.** It is a third data point. Do not substitute
it for `pipeline/baseline/build_baseline.sh`.

## Where it sits in the pipeline

```
HPRC assemblies ──┬─> PGGB monolithic ──────> baseline.gfa ─┐
                  │                                          ├─> Truvari
                  ├─> PGGB chunks ─> merge ──> merged.gfa ───┘   comparison
                  │
                  └─> dipcall / SVIM-asm ───> *.sv.vcf.gz ──> orthogonal check
```

## Tools

| Tool | Role | Alt allele style | Notes |
|---|---|---|---|
| **dipcall** | primary | sequence-resolved | Diploid-aware; hap1+hap2 jointly. Emits a confident-region BED. |
| **SVIM-asm** | alternative | symbolic (`<DEL>`, `<INV>`) | Detects inversions and tandem dups that dipcall does not type. |

Both are optional. If the toolchain is absent, both scripts print a
`SKIP` and exit 0 so the graph pipeline is never blocked.

## Install

```bash
conda install -c bioconda dipcall minimap2 k8 samtools htslib bcftools
conda install -c bioconda svim-asm          # optional alternative
```

## Usage

```bash
# All samples in work/manifests/hprc_selected.csv
bash pipeline/linear/run_dipcall.sh

# One sample, more threads, custom reference
REF=work/reference/GRCh38_chr21.fa THREADS=16 \
  bash pipeline/linear/run_dipcall.sh HG00673

# Alternative caller
bash pipeline/linear/run_svim_asm.sh HG00673
```

### Environment variables

| Variable | Default | Meaning |
|---|---|---|
| `REF` | `work/reference/GRCh38_chr21.fa` | Reference FASTA |
| `OUTDIR` | `results/linear/dipcall` | Output directory |
| `THREADS` | `min(nproc, 8)` | minimap2 threads |
| `MIN_SVLEN` | `50` | Minimum \|SVLEN\| for the SV subset |
| `DOWNLOAD_DIR` | `work/downloads` | Where assemblies were fetched |
| `MANIFEST` | `work/manifests/hprc_selected.csv` | Sample table |
| `MALE` / `PAR_BED` | `0` / unset | PAR handling for male samples on chrX/chrY |

## Haplotype convention

dipcall takes `hap1` then `hap2`. This repo passes **paternal as hap1 and
maternal as hap2**, matching the HPRC trio-phasing convention and the
`haplotype` numeric column in the official Release 2 index (1 = paternal,
2 = maternal in the manifest we derive). If you swap them the VCF is not
wrong, but the phase blocks will be inverted relative to the graph's
W-lines, and per-haplotype comparisons will silently mismatch.

## Outputs

```
results/linear/dipcall/
  HG00673.mak            # generated makefile (provenance: exact commands run)
  HG00673.dip.vcf.gz     # raw diploid calls, all sizes
  HG00673.dip.bed        # confident regions — use as Truvari --includebed
  HG00673.sv.vcf.gz      # PASS, biallelic, |SVLEN| >= MIN_SVLEN
results/linear/svim_asm/
  HG00673/variants.vcf   # raw SVIM-asm output
  HG00673.sv.vcf.gz      # sorted + indexed
results/logs/
  dipcall_HG00673.log
  svim_asm_HG00673.log
```

The `.dip.bed` confident-region file matters. Regions outside it are
where the assembly alignment was ambiguous, so calls there are unreliable
in *both* directions. Pass it to Truvari as `--includebed` when comparing,
or recall numbers will be penalized for regions no method can resolve.

## Comparing against graph calls

`pipeline/benchmark/benchmark_variants.sh` picks up
`results/linear/dipcall/*.sv.vcf.gz` automatically and runs a Truvari
comparison against the merged graph's deconstructed VCF.

Two representation caveats:

1. **SVIM-asm emits symbolic ALTs.** Truvari resolves `<DEL>`/`<INV>` only
   when `--reference` is supplied and the event is under `--max-resolve`
   (25 kbp). Always pass `--reference` and `--dup-to-ins` when comparing
   SVIM-asm output to sequence-resolved graph calls.
2. **dipcall emits sequence-resolved ALTs**, which match graph-derived
   calls much more cleanly. Prefer dipcall as the orthogonal reference.

## Interpreting disagreement

| Pattern | Likely cause |
|---|---|
| Graph finds it, dipcall does not | Real but assembly-alignment-ambiguous; or a graph construction artifact. Check whether it sits inside `.dip.bed`. |
| dipcall finds it, graph does not | PGGB segment length / identity thresholds filtered it, or it fell in a chunk boundary. |
| Both find it, different coordinates | Left-shift normalization. Re-run with `truvari refine --align mafft`. |
| Disagreement clusters at chunk boundaries | The reassembly stitch, not biology. This is the failure mode the project exists to detect. |

That last row is the one to watch. Chunk-boundary-correlated disagreement
is direct evidence about the research question, so
`benchmark_variants.sh` prints a hint pointing at `fn.vcf`/`fp.vcf` when
the verdict is `DIVERGENT`.

## Limitations

- Assembly-based calling inherits assembly errors; hifiasm misjoins become
  false SVs.
- Only variation representable relative to GRCh38 is visible. Sequence
  present in the samples but absent from the reference is invisible here,
  which is precisely the gap the pangenome graph is meant to close — so do
  not treat dipcall as ground truth, only as an independent witness.
- Small variants (< 50 bp) are excluded from the SV subset by default.
