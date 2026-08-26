# Parallel Pangenome Graph Construction (Linear Scaling SV)

> **Research question:** Can regional pangenome graphs be constructed independently in parallel and subsequently reassembled into a single graph while preserving haplotype paths, topology, sequence content, and variant representation relative to a conventionally constructed monolithic graph?

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
| Tests | ✅ 36 pass, 2 skip | GFA, chunking, merge, interval mapping, W-line |
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

## Repository Structure

```
pangenome-parallel/
|-- config/               # Pipeline configuration (YAML)
|-- dnanexus/applets/pggb_chunk/  # PGGB applet for parallel execution
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
|-- scripts/fetch_hprc_index.py  # Official HPRC Release 2 index
|-- scripts/download_hprc.py     # Download with gzip validation
|-- scripts/prepare_reference.sh  # GRCh38 chr21 from NCBI/DNAnexus
|-- scripts/setup_demo.py        # Synthetic end-to-end demo
|-- scripts/check_environment.sh  # Tiered environment checker
|-- tests/                       # 36 tests, pytest
|-- web/                         # Next.js + React + TypeScript
|-- docs/                        # Architecture, methods, DNAnexus
```

## License

MIT — see [LICENSE](LICENSE)

## Repository

GitHub: https://github.com/collaborativebioinformatics/linear-scaling-pangenome-sv
| HPRC index fetcher | ✅ | Official Release 2 index, match by assembly_name |
| Chromosome mapping | ✅ | minimap2-based interval mapping |
| Sequence extraction | ✅ | samtools faidx per-contig, no whole-genome RAM |
| DNAnexus staging | ✅ | Stage .fa.gz + .fai + .gzi + .md5 |
| DNAnexus applet | ✅ | pggb_chunk applet built via dx build |
| Pipeline orchestrator | ✅ | 9-checkpoint flow |
| Docker helpers | ✅ | PGGB, vg, odgi containers |
| Linear assembly pipeline | 📋 P1 | dipcall wrapper pending |
| Variant comparison | 📋 P1 | vg deconstruct + Truvari via Docker |