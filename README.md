# Parallel Pangenome Graph Construction (Linear Scaling SV)

> **Research question:** Can regional pangenome graphs be constructed independently in parallel and subsequently reassembled into a single graph while preserving haplotype paths, topology, sequence content, and variant representation relative to a conventionally constructed monolithic graph?

---

## Methods

<img src="slides/Screenshot 2026-08-27 at 8.53.22 AM.png" width="800"><br>
<img src="slides/Screenshot 2026-08-27 at 8.55.52 AM.png" width="800">
<img src="slides/Screenshot 2026-08-27 at 8.53.46 AM.png" width="800"><br>

Full editable deck: [`methods_slides.pptx`](methods_slides.pptx)

---

## Who Is Doing What?

| Person | Primary Area | What They Are Working On |
|--------|-------------|--------------------------|
| **Michael** | Pangenome graph construction | Lead PGGB/graph construction, graph method validation, biological graph troubleshooting |
| **Khoi** | Pipeline + linear/assembly workflow | Pipeline coding, linear-reference assembly-based SV calling, method validation |
| **Quang** | Pipeline architecture + graph merging | Pipeline architecture, chunk graph merging/reassembly algorithm, software engineering, testing |
| **Ali** | Pipeline + DNAnexus + integration + web | Pipeline coding, DNAnexus integration, benchmarking, end-to-end integration, web visualizer |
| **Alexander** | Documentation + presentation | Documentation, presentation, project story, explainability |

> **Important:** The pipeline is a shared technical responsibility. Michael, Khoi, Quang, and Ali collaboratively implement and troubleshoot the pipeline with overlapping ownership.

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
| Overlap-aware stitch | 🟡 NOT_IMPLEMENTED | Next algorithm milestone |
| Synthetic demo | ✅ | `make demo` generates full vertical slice |
| Tests | ✅ 55 pass, 2 skip | GFA, chunking, merge, interval mapping, W-line, provenance, parent-locus, orchestration |
| Web app | ✅ Initial | Next.js, loads JSON, dashboard |
| Environment checker | ✅ | Pipeline-tiered: REQUIRED / CONTAINER / WEB_OPTIONAL |
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
