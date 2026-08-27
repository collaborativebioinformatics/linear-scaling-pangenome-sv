
# Parallel Pangenome Graph Construction (Linear Scaling SV)

> **Research question:** Can regional pangenome graphs be constructed independently in parallel and subsequently reassembled into a single graph while preserving haplotype paths, topology, sequence content, and variant representation relative to a conventionally constructed monolithic graph?
## Quick Start
---
!!!! to be added soon :) !!!!
---

<img width="2347" height="1660" alt="logo" src="https://github.com/user-attachments/assets/e0a37556-7659-4a01-becb-7b6bf4b1e008" />

## Introduction

A pangenome graph represents the genetic variation across many individuals or haplotypes in a single data structure, rather than forcing every genome to be described as a list of differences from one linear reference sequence. Conventional linear comparisons use references like GRCh38, which collapse human diversity into a single representative path. With this approach, structural variants are often missed or misrepresented in regions of high diversity or repetitiveness. A pangenome graph can encode the divergent portions of a sample genome as alternate paths/nodes, inherently capturing structural variation such as large insertions, deletions, and rearrangements. 

The current standard approach to building these graphs, and the one this project builds on, is PGGB (the PanGenome Graph Builder). PGGB constructs a graph in three main stages: first, it performs all-pairs sequence alignment across the input haplotypes to identify homologous regions; second, it induces a variation graph from those pairwise alignments, collapsing shared sequence into common nodes, and lastly, it smooths and normalizes the graph through partial-order alignment so that regions that are locally linear are represented cleanly, while true structural variation remains visible as bubbles in the graph topology. The result is a single graph, plus a path for each input haplotype through that graph, that can be used for variant calling, alignment, and comparative analysis. The all-pairs alignment step compares every haplotype against every other haplotype, so its cost grows roughly quadratically (O(n²)) as more samples are added to the pangenome. Additional steps of the pipeline, alignment, induction, and smoothing, have to be processed across the full length of the input sequence, so cost also grows with the size of the genomic region being modeled. Thus, computation quickly becomes unwieldy and exceeds runtime and memory constraints of typical compute infrastructure. Our goal is to address the quadratic scaling in computation as sample size and length increases.

Our strategy is to break the monolithic construction problem into smaller, independent problems that can be solved in parallel and then reassembled. To accomplish this, we use minimap2 to identify haplotypes across our individual pangenomic samples. We then partition our samples into chunks that are aligned according to the haplotypes. Then, in parallel, we are able to build separate PGGB graphs for each chunk. We are able to re-stitch these graphs together across overlapping boundaries, assembling into a unified graph. 

To confirm that this parallel-and-stitch approach is actually a valid substitute for the conventional method, we validate the reassembled graph against a monolithic PGGB graph built the standard way over the same interval and haplotypes. Validation checks whether the two graphs agree on topology, whether each haplotype's path is preserved intact and contiguous across chunk boundaries, whether sequence content matches, and whether variant representation is equivalent. For that last point, we call variants from both graphs and compare them using truvari to measure precision, recall, and F1 between the merged-graph callset and the monolithic baseline. As a further, independent check, we generate a linear-reference SV callset for the same haplotypes against GRCh38 using assembly-based calling, giving us a second point of comparison outside the graph-based pipeline entirely. Together, these checks are meant to establish whether parallel construction can match conventional pangenome graph construction in accuracy while avoiding its computational bottleneck.



---

## Methods

For our pilot, we used four haplotype-resolved assemblies from HPRC Release 2 alongside GRCh38, restricted to a 1 Mb interval on chr21 (chr21:20,000,000–21,000,000) chosen as a representative smoke-test region. Because each haplotype assembly has its own independent coordinate system, we first used minimap2 to align each haplotype to GRCh38 and project the target interval's boundaries onto each haplotype's own sequence, filtering low-confidence mappings below a minimum mapping quality of 20 and a minimum query coverage of 0.90; samtools-based interval extraction and a custom multi-FASTA builder then produced the corresponding per-haplotype sequence for the interval. From this shared input, we built two graphs for comparison: a monolithic baseline, produced by running PGGB v0.6.0 once over the full interval across all haplotypes in a Dockerized environment, and a parallel graph, produced by partitioning the same interval into overlapping 400 kb windows (50 kb pairwise overlap) with a custom chunking script, extracting each haplotype's corresponding sub-sequence per chunk under a parent-locus constraint that rejects chunks spanning contig boundaries or ambiguous strand orientation, and running PGGB independently and in parallel on each chunk with identical parameters. Both the local synthetic pipeline and full-scale runs on real HPRC data are orchestrated the same way: a Makefile drives local execution, while a matching set of DNAnexus applets and shell scripts stage inputs and launch the equivalent chunked and baseline jobs on cloud compute for real-data runs.

The independently built chunk graphs are represented and manipulated using a custom GFA data model and parser we implemented in-house, supporting the full GFA1 structure (headers, segments, links, paths, and the official W-line walk format) with round-trip parsing, serialization, and file I/O. Reassembly currently proceeds in two stages: a namespace-safe disjoint union that concatenates chunk graphs without merging shared content, used as a diagnostic baseline, and an overlap-aware stitching algorithm — the project's core methodological contribution, currently under active development — that identifies each haplotype's common subpath within a chunk-pair's overlap region and welds matching segment IDs across the boundary to produce a single continuous graph, with a boundary report tracking per-haplotype stitch success and base-level identity across each seam. For validation, we run vg deconstruct on both the baseline and reassembled graphs to call variants from each, then compare the two callsets with truvari to obtain precision, recall, and F1; independently, we generate an assembly-based linear-reference SV callset for the same haplotypes against GRCh38 using dipcall, giving us a comparison point outside the graph-based pipeline entirely. The full pipeline is covered by a pytest suite (55 tests spanning the GFA model, chunking, merging, interval mapping, and pipeline orchestration), and results are exported as compact JSON and served through a Next.js/React/TypeScript web dashboard for inspection.

<img src="slides/Screenshot 2026-08-27 at 8.53.22 AM.png" width="800"><br>
<img src="slides/Screenshot 2026-08-27 at 8.55.52 AM.png" width="800">
<img src="slides/Screenshot 2026-08-27 at 8.53.46 AM.png" width="800"><br>

Full editable deck: [`methods_slides.pptx`](methods_slides.pptx)

---

## Who Is Doing What?

**Full task breakdown with commands, file paths, and acceptance criteria:**
👉 **[docs/TASK_BREAKDOWN.md](docs/TASK_BREAKDOWN.md)**

| Person | Primary Area | What They Are Working On |
|--------|-------------|--------------------------|
| **Quang** | Pipeline architecture + graph merging | 🔴 **Q1 (HARD BLOCKER):** Overlap-aware stitch integration. Fix `setup_demo.py` config reading, `_load_chunks` path, and `chunk_rows` parameter. Merge `quang-overlap-aware-merge` branch to main. |
| **Michael** | Pangenome graph construction | M1: Run real PGGB applet on DNAnexus (micro-test). M2/M3: Build baseline + chunk graphs on DNAnexus. M4: Full 1 Mb smoke test. |
| **Khoi** | Pipeline + linear/assembly workflow | K1: Dipcall wrapper validation. K2: SVIM-asm wrapper. K3: Truvari variant comparison. K4: Method documentation. |
| **Ali** | Pipeline + DNAnexus + integration + web | ✅ A1-A5 complete. 🔵 **A6-A10 (current):** Web visualization — data export pipeline, benchmark display, chunk visualization, compare mode, Vercel deploy. |
| **Alexander** | Documentation + presentation | D1: Update architecture doc. D2: DNAnexus operations guide. D3: Web visualization guide. D4: Final presentation. |
| **Lex** | PGGB benchmarking | L1: Local PGGB cost benchmarks. L2: Graph cost model. L3: Document findings. L4: Merge into pipeline. |

> **Important:** The pipeline is a shared technical responsibility. Everyone collaboratively implements and troubleshoots.

---

## Where Should I Work?

```
Want to work on PGGB / graph construction?
  pipeline/baseline/ and pipeline/parallel/  (Michael + Khoi)

Want to work on graph stitching?
  pipeline/merge/  (Quang + Michael + Khoi + Ali)

Want to work on DNAnexus / cloud execution?
  dnanexus/  (Ali + Quang + Khoi)

Want to work on benchmarking / validation?
  pipeline/benchmark/  (Ali + Khoi + Michael)

Want to work on the website / visualization?
  web/  (Ali + anyone)

Want to work on writing / presentation?
  docs/  (Alexander)
```

---

## Pipeline Flow

```
HPRC Release 2 (4 haplotypes + GRCh38)
  -> minimap2 interval mapping
  -> monolithic PGGB (baseline)  OR  chunk PGGB (parallel)
  -> merge (NOT_IMPLEMENTED - disjoint union only)
  -> validation & benchmark
  -> JSON -> web
```

**Merge status:** NOT_IMPLEMENTED. Current merged graph is disjoint union only (diagnostic). Overlap-aware stitching is the next algorithm milestone.

---

## Current State

| Component | Status | Notes |
|-----------|--------|-------|
| GFA data model | ✅ | Header, Segment, Link, Path, Walk (official W-line, no step_count) |
| GFA parser | ✅ | Parse, dump, roundtrip, file I/O |
| Chunking | ✅ | Overlapping chunks, configurable pairwise overlap |
| Merge (disjoint union) | ✅ | Namespace-safe concatenation, diagnostic only |
| Overlap-aware stitch | 🟡 NOT_IMPLEMENTED | Algorithm exists on quang branch, needs integration (Q1) |
| Synthetic demo | ✅ | `make demo` generates full vertical slice |
| Tests | ✅ 105 pass, 2 skip | GFA, chunking, merge, interval mapping, W-line, provenance, graph stats |
| Web explorer | ✅ Interactive | Cytoscape.js graph, sample/haplotype selector, inspector, 4 tabs |
| Web JSON guard | ✅ | `guard_no_genomic()` blocks GFA/FASTA/VCF |
| DNAnexus applets | ✅ Built | `pggb_chunk` + `pggb_baseline`, instance `mem3_ssd1_v2_x16` |
| DNAnexus dry-run | ✅ | 3 chunks validated, FATAL checks verified |
| Graph statistics | ✅ | Full topology, comparison, TSV export |
| Environment checker | ✅ | Pipeline-tiered: REQUIRED / CONTAINER / WEB_OPTIONAL |
| Stitch boundary validation | 🟡 NOT_RUN | Requires Q1 first |
| Equivalence validation | 🟡 NOT_RUN | Requires Q1 + M1-M4 first |
| Variant comparison | 🟡 NOT_RUN | Requires K1-K3 first |
---

## Quick Start

### Local (synthetic demo)
```bash
make check      # Check environment
make demo       # Generate synthetic data, build GFAs, merge, export JSON
make test       # Run 36+ tests
make web        # Start Next.js dev server at http://localhost:3000
```

### DNAnexus (real HPRC data)
```bash
# Inside DNAnexus Cloud Workstation:
cd ~/pangenome-parallel
git pull origin main
bash dnanexus/setup_workstation.sh   # One-time setup
bash dnanexus/run_pipeline.sh --upload  # Full pipeline
```

## Configuration

| Parameter | File | Default | Description |
|-----------|------|---------|-------------|
| target.chromosome | config/pipeline.yaml | chr21 | Target chromosome |
| target.start/end | config/pipeline.yaml | 20000000-21000000 | Smoke interval (0-based half-open) |
| mapping.min_query_coverage | config/pipeline.yaml | 0.90 | Minimum query coverage |
| mapping.min_mapq | config/pipeline.yaml | 20 | Minimum MAPQ |
| parallel.chunk_size_bp | config/pipeline.yaml | 400000 | Chunk size in bp |
| parallel.overlap_bp | config/pipeline.yaml | 50000 | Pairwise overlap between chunks |
| pggb.threads | config/pipeline.yaml | 8 | PGGB thread count |

## What To Do Next (Pre-Submission)

All 19 P0 audit items are fixed. Below is the exact pre-submission task list broken down per person. **Do not launch PGGB on real data until the stitch is implemented.**

### Quang — Merge Algorithm (P0, highest priority)

Stitch is the single critical unimplemented component.

| # | Task | Area | Instructions |
|---|------|------|-------------|
| Q1 | **Overlap-aware stitch** | `pipeline/merge/merge_graphs.py` | Replace `NOT_IMPLEMENTED` path with real `overlap_aware_stitch()`. Algorithm: (1) load two adjacent chunk GFAs, (2) find common subpath in the overlap region per haplotype, (3) weld by merging shared segment IDs, (4) write stitched output to `results/merge/merged.gfa`. |
| Q2 | **Stitch tests** | `tests/test_merge.py` | Test stitch on small overlapping GFAs: path continuity, sequence identity across boundary, no orphan nodes. |
| Q3 | **Boundary report** | `pipeline/merge/validate_merge.py` | Per-haplotype stitch success/failure, base-level identity across boundary. |

### Michael — PGGB + Graph Validation (P0)

| # | Task | Area | Instructions |
|---|------|------|-------------|
| M1 | **PGGB version validation** | `config/pipeline.yaml` | Verify the pinned image digest matches actual `v0.6.0` tag on ghcr.io. Update if mismatch. |
| M2 | **Chunk FASTA review** | `pipeline/parallel/build_all_chunks.py` | Review parent-locus constraint (lines 180-202): contig FATAL, strand warning, margin calc. |
| M3 | **Baseline graph review** | `pipeline/baseline/build_baseline.sh` | Confirm `-w 50000` used by both baseline and chunks. Verify command matches manual pggb invocation. |
| M4 | **Test data creation** | `tests/` | Create a small (5 Kb) test FASTA with 2-3 haplotypes for stitch testing without real HPRC data. |

### Khoi — Linear Pipeline + Validation (P0/P1)

| # | Task | Area | Instructions |
|---|------|------|-------------|
| K1 | **Linear dipcall wrapper** | `pipeline/linear/run_dipcall.sh` | Download/verify dipcall, run on 4 HPRC haplotypes vs GRCh38 chr21, produce VCF. Needed for benchmark. |
| K2 | **Variant comparison** | `pipeline/benchmark/benchmark_variants.sh` | `vg deconstruct` on baseline + merged GFAs, `truvari` compare baseline vs merged VCF. Record precision/recall/F1. |
| K3 | **Benchmark audit** | `pipeline/benchmark/benchmark_variants.sh` | Fix existing stub. Wire results into `results/benchmark/`. |
| K4 | **Real data smoke test** | — | Once stitch ready (Q1), run full pipeline on chr21:20000000-21000000. Report baseline/chunks/stitch/equivalence status. |

### Ali — DNAnexus + Integration + Web (P0/P1)

| # | Task | Area | Instructions |
|---|------|------|-------------|
| A1 | **Baseline applet test** | `dnanexus/applets/pggb_baseline/` | Build: `dx build --brief`. Verify same instance type as chunk (`mem3_ssd1_v2_x16`). |
| A2 | **Orchestration dry run** | `dnanexus/run_parallel_chunks.sh` | Verify expected=local=uploaded=submitted=downloaded. Simulate missing FASTA, verify FATAL. |
| A3 | **Timing dashboard** | `dnanexus/run_parallel_chunks.sh` | Review `graph_parallel_wall_seconds` (max stop - min start) and `sum_worker_seconds`. |
| A4 | **Web stitch status** | `web/app/page.tsx` | Show `stitch: NOT_IMPLEMENTED` badge. Label merged graph as "diagnostic only". |
| A5 | **Web JSON guard** | `scripts/sync_web_results.py` | Prevent accidental GFA copy into `web/public/data/`. Only compact JSON allowed. |

### Alexander — Documentation + Presentation (P0)

| # | Task | Area | Instructions |
|---|------|------|-------------|
| D1 | **Architecture diagram** | `docs/ARCHITECTURE.md` | Draw pipeline architecture (ASCII/Mermaid). Highlight merge as NOT_IMPLEMENTED. |
| D2 | **Methods summary** | `docs/` | 1-2 pages: parallelization strategy, parent-locus constraint, provenance hashing, benchmark design. |
| D3 | **Presentation slides** | `docs/presentation/` | Research question, architecture, status, timeline, next steps. |
| D4 | **README review** | `README.md` | Verify "Who Is Doing What" table. Read for clarity. |

### Shared — Critical Path

| # | Task | Who | Instructions |
|---|------|-----|-------------|
| S1 | **Integration test** | Quang + Michael + Ali | After Q1, run `make demo`. Verify baseline builds, chunks build, stitch merges, paths match. |
| S2 | **Final validation** | Everyone | `make check && make test` from clean checkout. All 55+ tests pass. CI on the merge commit. |
| S3 | **Tag release** | Ali | Squash-merge working branch to main. Tag `v0.1.0-pre`. Do NOT launch real PGGB. |

### Timeline

```
Now ___ Q1 (stitch: 1-2d) ___ S1 (integration: 1d) ___ S2 (validation: 1d) ___ Submission
        M1-M4, A1-A5 (parallel)   D1-D4 (parallel)        K1-K4 (parallel)
```

### Stitch Algorithm Reference (for Q1)

1. **Input:** Two adjacent chunk GFAs with known overlap region from `chunk_manifest.tsv`.
2. **Find common subpaths:** Extract each haplotype's path segment covering the overlap.
3. **Create weld mapping:** Map segment IDs from chunk_B's overlap to chunk_A's by sequence identity.
4. **Merge nodes:** Replace chunk_B IDs with chunk_A IDs in links/paths. Drop duplicate segments.
5. **Concatenate remaining:** Keep both graphs' content outside the overlap.
6. **Validate:** For every haplotype, verify path continuity and sequence identity to monolithic baseline.

Start with exact-match welding before attempting alignment-based welding for near-matches.

### Submission Checklist

```
☐ Q1: Stitch implemented and tested
☐ Q2: Stitch tests pass
☐ M1: PGGB version verified
☐ K1: Dipcall wrapper complete
☐ K2: Variant comparison working
☐ A1-A5: DNAnexus + web updates done
☐ D1-D4: Documentation complete
☐ S1: Full integration test passes
☐ S2: All 55+ tests pass on CI
☐ S3: Tagged v0.1.0-pre
```

## Repository Structure

```
pangenome-parallel/
|-- config/               # Pipeline configuration (YAML)
|-- dnanexus/applets/pggb_chunk/  # PGGB applet for parallel chunks
|-- dnanexus/applets/pggb_baseline/  # PGGB applet for baseline
|-- dnanexus/run_pipeline.sh      # Master orchestrator
|-- dnanexus/run_parallel_chunks.sh  # Launch chunk jobs on DNAnexus
|-- dnanexus/stage_inputs.sh      # Stage HPRC assemblies from storage
|-- dnanexus/setup_workstation.sh  # Bootstrap Cloud Workstation
|-- pipeline/merge/gfa.py         # GFA data model (core library)
|-- pipeline/merge/merge_graphs.py  # Merge strategies
|-- pipeline/parallel/make_chunks.py  # Overlapping chunk creation
|-- pipeline/parallel/build_all_chunks.py  # Chunk FASTA builder
|-- pipeline/baseline/build_baseline.sh  # Monolithic PGGB via Docker
|-- pipeline/prepare/faidx_utils.py  # samtools interval extraction
|-- pipeline/prepare/map_chromosome.py  # minimap2 interval mapping
|-- pipeline/prepare/prepare_sequences.py  # Multi-FASTA builder
|-- pipeline/benchmark/           # Stats, path compare, reports
|-- scripts/run_pggb.py           # Canonical PGGB runner (reads config)
|-- scripts/fetch_hprc_index.py  # Official HPRC Release 2 index
|-- scripts/download_hprc.py     # Download with gzip validation
|-- scripts/prepare_reference.sh  # GRCh38 chr21 from DNAnexus/NCBI
|-- scripts/setup_demo.py        # Synthetic end-to-end demo
|-- scripts/check_environment.sh  # Tiered environment checker
|-- tests/                       # 55 tests, pytest
|-- web/                         # Next.js + React + TypeScript
|-- docs/                        # Architecture, methods, DNAnexus
```

## License

MIT — see [LICENSE](LICENSE)

## Repository

GitHub: https://github.com/collaborativebioinformatics/linear-scaling-pangenome-sv
