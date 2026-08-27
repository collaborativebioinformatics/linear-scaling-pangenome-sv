# Team Task Breakdown — Final Sprint to Hackathon Demo

> Ali is now focused exclusively on web/visualization (A6-A10).
> This document details every remaining task with precise file paths, commands, and acceptance criteria.

---

## Critical Dependency Chain

```
Quang Q1 (stitch) -> Michael M1 (PGGB baseline) -> Khoi K1 (dipcall)
     |                        |
     v                        v
Lex L1 (benchmarking)    Everyone S1 (integration)
     |
     v
Ali web (data pipeline) -> Alexander D1 (documentation)
```

**Quang Q1 is the hard blocker.** Until the overlap-aware stitch is integrated, the merged graph is a disjoint union and equivalence cannot be validated.

---

## Q1 — Quang: Overlap-Aware Stitch (HARD BLOCKER)

### Current State
- Algorithm exists on `origin/quang-overlap-aware-merge` branch
- `pipeline/merge/merge_graphs.py` has `overlap_aware_stitch()` function
- `pipeline/merge/gfa.py` on main has correct W-line parser + all 11 helpers
- **Two bugs blocking integration** (from MERGE_NOTES.md)

### Bug 1: `setup_demo.py` ignores `config.merge.strategy`
**File:** `scripts/setup_demo.py` line ~117
**Fix:** Read `config.merge.strategy` instead of hardcoding `diagnostic_disjoint_union()`

### Bug 2: `_load_chunks` hardcodes path
**File:** `pipeline/merge/merge_graphs.py`
**Fix:** Make `chunk_dir` a parameter defaulting to `work/chunks`

### Bug 3: `overlap_aware_stitch` requires `chunk_rows`
**File:** `pipeline/merge/merge_graphs.py`
**Fix:** Auto-discover from `chunk_manifest.tsv` when not provided

### Acceptance Criteria
```bash
make demo
# Expected: "Strategy: overlap_aware" in output
# Expected: 5 components (not 15), EQUIVALENT verdict
```

### Merge Instructions
```bash
git fetch origin quang-overlap-aware-merge
git checkout -b quang-merge-fix origin/quang-overlap-aware-merge
# Fix bugs 1-3 above
git add -A && git commit -m "fix: merge strategy, chunk paths, optional chunk_rows"
git push origin quang-merge-fix
# PR to main
---

## M1-M4 — Michael: PGGB Graph Construction

### M1: Verify PGGB v0.6.0 on DNAnexus
**Why:** Applets are built (A1) but never actually executed.
**Command:**
```bash
cd ~/pangenome-parallel
git pull origin main
dx upload tests/fixtures/tiny_reference.fa --destination /data/test/
PGGB_CONFIG=$(python3 scripts/gen_pggb_config.py)
dx run applet-JB88PZ00ZQv43jJPQ11v2pVp \
  -i fasta=/data/test/tiny_reference.fa \
  -i pggb_config_json="$PGGB_CONFIG" \
  --instance-type mem3_ssd1_v2_x16 \
  --destination /results/microtest/ \
  --name "Micro PGGB Test" --brief
```
**Acceptance:** GFA > 0 bytes, metadata has config_sha256, log shows PGGB v0.6.0.

### M2: Build Baseline Graph on DNAnexus
**File:** `pipeline/baseline/build_baseline.sh`
**Acceptance:** baseline.gfa, baseline.log, baseline_metadata.json all present.

### M3: Build Chunk Graphs on DNAnexus
**File:** `dnanexus/run_parallel_chunks.sh`
**Acceptance:** All chunks succeed, expected == uploaded == submitted == downloaded.

### M4: Run Full 1 Mb Smoke Test
**Command:** `make check && dnanexus/run_pipeline.sh`
**Acceptance:** Baseline builds, chunks build, stitch merges, paths continuous.

---

## K1-K4 — Khoi: Linear/Assembly Pipeline

### K1: Dipcall Wrapper
**File:** `pipeline/linear/run_dipcall.sh`
**Purpose:** Run dipcall to get truth VCF for chr21 region.
**Acceptance:** VCF output with variants against GRCh38 chr21:20000000-21000000.

### K2: SVIM-asm Wrapper
**File:** `pipeline/linear/run_svim_asm.sh` (now executable)
**Acceptance:** VCF output.

### K3: Variant Comparison
**File:** `pipeline/benchmark/benchmark_variants.sh`
**Acceptance:** Truvari summary JSON with precision/recall.

### K4: Method Validation
**Deliverable:** `docs/methods/linear-pipeline.md`

---

## L1-L4 — Lex: PGGB Benchmarking

### L1: Run Local PGGB Cost Benchmark
**File:** `lex_testing/pggb_run.sh`
**Acceptance:** Benchmark CSV with wall time, peak memory, CPU.

### L2: Graph Cost Model
**File:** `lex_testing/benchmark.sh`
**Acceptance:** Plot showing scaling behavior.

### L3: Document Findings
**Deliverable:** `docs/benchmarking/pggb-scaling.md`

### L4: Merge Benchmarking into Pipeline
**Acceptance:** lex_testing results in results/benchmark/ for web display.

---

## D1-D4 — Alexander: Documentation & Presentation

### D1: Architecture Document
**File:** `docs/ARCHITECTURE.md` (update existing)

### D2: DNAnexus Operations Guide
**File:** `docs/DNANEXUS.md` (update existing)

### D3: Web Visualization Guide
**File:** `docs/web-visualization.md` (new)

### D4: Final Presentation
**Slides:** `presentation/` directory

---

## S1-S3 — Everyone: Integration & Submission

### S1: Full Integration Test
**When:** After Q1, M1-M4, K1-K2 are all done.
**Command:** `make check && make test && make demo`
**Acceptance:** All checks pass, 105+ tests pass, EQUIVALENT verdict.

### S2: CI Pipeline
**File:** `.github/workflows/ci.yml`

### S3: Tag Release
```bash
git tag -a v0.1.0-pre -m "Hackathon pre-release"
git push origin v0.1.0-pre
```

---

## A6-A10 — Ali: Web Visualization (My Focus)

### Current State
- ✅ Cytoscape.js interactive graph with zoom/pan/search
- ✅ Sample/haplotype selector sidebar
- ✅ Inspector panel for node/edge details
- ✅ Tab layout: Explore, Dashboard, Chunks, Compare
- ✅ Synthetic demo graph generation
- ✅ Pipeline status badges
- ✅ TypeScript types
- ✅ Web typecheck + build passing

### A6: Real Data Pipeline (new file)
**File:** `pipeline/export/build_web_dataset.py`
**Purpose:** Convert GFA to compact JSON for browser.
**Design:** One JSON per sample/haplotype, max 1500 nodes/4000 edges.

### A7: Benchmark Data Display
**Purpose:** Show real benchmark metrics when available.

### A8: Chunk Visualization
**Purpose:** Chunk overlap diagram with coordinate ranges.

### A9: Compare Mode
**Purpose:** Side-by-side haplotype comparison.

### A10: Deploy to Vercel
**Command:** `cd web && npx vercel --prod`

---

## Quick Reference: Key Commands

```bash
# Run everything locally (synthetic)
make check && make test && make demo
cd web && npm run dev

# Build applets on DNAnexus
cd dnanexus/applets/pggb_chunk && dx build --brief --destination /applets/pggb_chunk/
cd ../pggb_baseline && dx build --brief --destination /applets/pggb_baseline/

# Run pipeline
bash dnanexus/run_pipeline.sh

# Deploy web
cd web && npm run build && npx vercel --prod
```
```
