# Parallel Pangenome Graph Construction

---

## Who Is Doing What?

| Person | Primary Area | What They Are Working On |
|--------|-------------|--------------------------|
| **Michael** | Pangenome graph construction | Lead PGGB/graph construction, help write graph code, verify graph-building methods and parameters match accepted/current pangenome practices, troubleshoot biological graph issues |
| **Khoi** | Pipeline + linear/assembly workflow | Actively code the pipeline with the team, implement/validate the linear-reference assembly-based workflow, investigate papers/methods while implementing, help troubleshoot graph construction |
| **Quang** | Pipeline architecture + graph merging | Build pipeline architecture with the team, develop the chunk graph merging/reassembly algorithm, software engineering, testing and integration |
| **Ali** | Pipeline + DNAnexus + integration + web | Work directly alongside Khoi, Michael, and Quang on pipeline code; DNAnexus integration; benchmarking; end-to-end integration; troubleshooting; build the web visualizer and Vercel deployment |
| **Alexander** | Documentation + presentation | Lead documentation, presentation, project story, and explainability; help the team reason about overall system design |

> **Important:** The pipeline is a shared technical responsibility. Ali is not building the pipeline alone. Michael, Khoi, Quang, and Ali should actively collaborate on implementation and troubleshooting, with overlapping ownership so no critical technical component has only one person who understands it.

---

## Your Code: What to Review & Improve

Each person has files marked with `TODO (YourName)` throughout the codebase. Here is your specific review list:

### 👨‍🔬 Michael — Graph Construction & PGGB

| File | What Needs You |
|------|---------------|
| `pipeline/baseline/build_baseline.sh` | Validate PGGB parameters for chr21 HPRC data (line 30-45: -p 90, -s 5000, -k 29). Add container fallback. |
| `pipeline/parallel/build_chunk.sh` | Same PGGB parameter validation at chunk scale. Add chunk-specific logging. |
| `pipeline/merge/gfa.py` | Review `infer_data_mode()` heuristic — is it correct for real vs synthetic? Add more graph validation methods. |
| `pipeline/prepare/map_chromosome.py` | This is critical: implement the contig-to-chromosome mapping using minimap2 or HPRC metadata. Never assume coordinates match. |
| `pipeline/prepare/extract_regions.py` | Implement sequence extraction using mapped coordinates. |
| `tests/test_gfa_parser.py` | Add tests for edge cases: empty sequences, negative strands, cyclic graphs, large node IDs. |

### 👨‍💻 Khoi — Pipeline & Linear Assembly

| File | What Needs You |
|------|---------------|
| `pipeline/linear/run_dipcall.sh` | Implement actual dipcall wrapper for HG00673 and HG00733. Add SVIM-asm as alternative. |
| `pipeline/linear/README.md` | Document the assembly-based SV workflow in more detail. |
| `scripts/fetch_hprc_index.py` | Verify HPRC Release 2 index URL still works; add retry logic. |
| `pipeline/benchmark/benchmark_variants.sh` | Implement vg deconstruct + Truvari comparison when tools are available. |
| `pipeline/benchmark/graph_stats.py` | Add component analysis, branching node count, N50-style metrics. |
| `pipeline/prepare/prepare_sequences.py` | Add actual sequence extraction logic once assemblies are downloaded. |

### 👨‍🎓 Quang — Graph Merging & Pipeline Architecture

| File | What Needs You |
|------|---------------|
| `pipeline/merge/merge_graphs.py` | **This is the experimental core.** Implement the real overlap-aware stitch in `overlap_aware_stitch()`. Current version is a disjoint union placeholder. |
| `pipeline/merge/gfa.py` | Add `get_path_sequence()`, `get_walk_sequence()`, and other path traversal methods needed for merge validation. |
| `pipeline/merge/paths.py` | Implement path traversal across chunk boundaries. |
| `pipeline/merge/validate_merge.py` | Add more validation: check all samples preserved, reference continuity, no orphaned nodes. |
| `pipeline/parallel/make_chunks.py` | Review overlap coordinate math — is `overlap_left`/`overlap_right` computed correctly for edge chunks? |
| `pipeline/benchmark/compare_paths.py` | Add sequence-level comparison of spelled paths. |
| `tests/test_gfa_parser.py` | Add test for merge operations: disjoint union on real data, path preservation. |

### 👨‍🔧 Ali — Integration, Web, DNAnexus & End-to-End

| File | What Needs You |
|------|---------------|
| `web/` | Deploy to Vercel. Ensure `/api/data` route works. Add graph viewer, boundary diff view, pipeline status flow. |
| `dnanexus/` | Implement actual `dx` commands in setup scripts. Add applets for PGGB execution. |
| `scripts/sync_web_results.py` | Wire up to copy real pipeline results to `web/public/data/`. |
| `pipeline/benchmark/build_report.py` | Add runtime comparison when PGGB run data is available. |
| `pipeline/export/gfa_to_json.py` | Add graph windowing for large GFAs so the web viewer doesn't crash. |
| `pipeline/export/find_bubbles.py` | Improve bubble detection and integrate with web display. |
| `Makefile` | Add Vercel deploy target (`make deploy`). Polish commands. |

### 👨‍🏫 Alexander — Documentation & Presentation

| File | What Needs You |
|------|---------------|
| `README.md` | Keep this updated as the project evolves. Add screenshots when the web app has real data. |
| `docs/ARCHITECTURE.md` | Verify architecture diagram matches actual implementation. Add component descriptions. |
| `docs/METHODS.md` | Write a clear methods section suitable for judges who may not be graph experts. |
| `docs/DEMO.md` | Create a step-by-step demo walkthrough for the hackathon presentation. |
| `pilot_pipeline.pptx` | Update with current results and architecture diagrams. |

---

## Where Should I Work?

```
Want to work on PGGB?
→ pipeline/baseline/ and pipeline/parallel/
   Michael + Khoi + anyone helping

Want to work on linear assembly-based calling?
→ pipeline/linear/
   Khoi + Michael

Want to work on graph stitching?
→ pipeline/merge/
   Quang + Khoi + Michael + Ali

Want to work on DNAnexus?
→ dnanexus/
   Ali + Quang + Khoi

Want to work on benchmarking?
→ pipeline/benchmark/
   Ali + Khoi + Michael

Want to work on the website?
→ web/
   Ali + anyone available

Want to work on writing/presentation?
→ docs/
   Alexander
```

---

## Today's Priority

```
P0 — BUILD THE PIPELINE

HPRC data → chr21 preparation → baseline PGGB graph → parallel chunks
→ chunk PGGB graphs → graph reassembly → validation → metrics/JSON → web visualizer
```

The goal today is to get the **linear workflow, graph-building workflow, parallel/merging algorithm, and web interface moving simultaneously**, connecting each component as soon as it produces usable output.
---

## What's Built So Far

### Repository Structure

```
pangenome-parallel/
├── README.md                  ← You are here
├── CONTRIBUTING.md            ← How to collaborate
├── Makefile                   ← make demo, make test, make web
├── .gitignore                 ← Ignores large files
├── config/                    ← Pipeline configuration (YAML)
├── pipeline/
│   ├── merge/gfa.py           ← GFA data model (core library)
│   ├── merge/merge_graphs.py  ← Merge strategies
│   ├── parallel/make_chunks.py← Chunk creation
│   ├── benchmark/             ← Stats, path compare, reports
│   ├── export/                ← GFA→JSON, bubble detection
│   └── prepare/               ← Sequence preparation
├── scripts/                   ← Setup, fetch, download, sync
├── tests/                     ← 17 passing tests
├── web/                       ← Next.js + React + TypeScript
├── docs/                      ← Architecture, methods, DNAnexus
├── dnanexus/                  ← Cloud deployment
├── results/                   ← Outputs (gitignored)
└── work/                      ← Working files (gitignored)
```

### Current State

| Component | Status | Notes |
|-----------|--------|-------|
| GFA data model | ✅ Working | Header, Segment, Link, Path, Walk, GfaGraph |
| GFA parser | ✅ Working | Parse, dump, roundtrip, file I/O |
| Chunking | ✅ Working | Overlapping chunk manifest creation |
| Merge (disjoint union) | ✅ Working | Namespace-safe concatenation |
| Synthetic demo | ✅ Working | `make demo` generates full vertical slice |
| Tests | ✅ 17 passing | GFA parsing, chunking, merge, export |
| Web app | ✅ Initial | Next.js, loads JSON, shows dashboard |
| Environment checker | ✅ Working | `bash scripts/check_environment.sh` |
| HPRC index fetcher | ✅ Implemented | Fetches from official HPRC index |
| Overlap-aware stitch | 🟡 Skeleton | Real implementation pending PGGB output |
| DNAnexus integration | 🟡 Started | Docs created, helpers pending |
| Linear assembly pipeline | 📋 Planned | dipcall wrapper pending |
| Benchmarking | 🟡 Framework | Stats, path compare, report built |
---

## Quick Start

```bash
# Check environment
bash scripts/check_environment.sh

# Run the synthetic demo (no external tools required)
make demo

# Run tests
make test

# Start the web app
make web
# OR: cd web && npm run dev
```

Open [http://localhost:3000](http://localhost:3000) to see the explorer.

## Data Pipeline Commands

```bash
make fetch-index     # Fetch HPRC Release 2 assembly index
make download        # Download selected HPRC assemblies
make prepare         # Prepare chr21 sequences
make baseline        # Build monolithic PGGB graph
make chunks          # Create chunk manifest + build chunk graphs
make merge           # Merge chunk graphs
make benchmark       # Run validation and benchmarking
make smoke           # Full smoke-test pipeline (HPRC chr21:20m-21m)
make chr21           # Full chr21 pipeline (requires DNAnexus)
make clean           # Remove results/ and work/
```

## Success Criteria

### P0 — today
- [x] Repository builds/tests locally
- [x] Synthetic demo works end to end
- [x] HPRC sample manifest can be generated from Release 2 index
- [x] DNAnexus workstation instructions documented
- [ ] Real HPRC chr21 smoke-test dataset prepared
- [ ] Monolithic PGGB graph produced
- [ ] Multiple chunk graphs produced
- [x] Merge code can consume chunk GFAs
- [x] Expected haplotype paths survive merge
- [x] Baseline and merged graph stats generated
- [x] Web application loads exported JSON

### P1
- [ ] Overlap-aware merge works across all smoke-test boundaries
- [ ] Variant comparison runs
- [ ] Runtime/memory benchmark runs
- [ ] Web app highlights boundary disagreements
- [ ] Real results deployed to Vercel

### P2
- [ ] Run complete chromosome 21
- [ ] Run assembly-based `dipcall`
- [ ] Compare graph-derived and assembly-derived calls
- [ ] Formal DNAnexus applets/workflows
- [ ] Test additional haplotypes/chromosomes

---

## License

MIT — see [LICENSE](LICENSE)

## Team

BCM SV Hackathon — 2026
> **Research question:** Can regional pangenome graphs be constructed independently in parallel and subsequently reassembled into a single graph while preserving haplotype paths, topology, sequence content, and variant representation relative to a conventionally constructed monolithic graph?