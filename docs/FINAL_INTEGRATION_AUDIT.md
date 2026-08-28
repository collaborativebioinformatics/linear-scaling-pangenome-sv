# Final Integration Audit

> Audit date: 2026-08-28
> HEAD at audit: `cf8055c` (before) → see final commit (after)

## Verification Results

| Check | Result |
|-------|--------|
| pytest | **148 passed, 2 skipped** |
| py_compile (pipeline, scripts, tests) | PASS |
| bash -n (all .sh) | PASS |
| dxapp.json (json.tool) | PASS |
| web typecheck (tsc --noEmit) | PASS |
| web build (next build) | PASS |
| make demo (synthetic end-to-end) | PASS — `overlap_aware`, EQUIVALENT, 2/2 boundaries |

## Component Status Table

Legend: **Implemented** = code exists and wired; **Static** = source/type checks;
**Synthetic** = synthetic demo validated; **DNAnexus** = real execution;
**Real HPRC** = real HPRC data validated.

| Component | Implemented | Static | Synthetic | DNAnexus | Real HPRC | Status |
|-----------|-------------|--------|-----------|----------|-----------|--------|
| Reference preparation | PASS | PASS | NOT_RUN | NOT_VERIFIED | NOT_RUN | PARTIAL |
| HPRC staging | PASS | PASS | NOT_RUN | NOT_VERIFIED | NOT_RUN | PARTIAL |
| Mapping | PASS | PASS | NOT_RUN | NOT_VERIFIED | NOT_RUN | PARTIAL |
| Chunk building | PASS | PASS | PASS | NOT_VERIFIED | NOT_RUN | PARTIAL |
| Baseline PGGB | PASS | PASS | PASS | NOT_VERIFIED | NOT_RUN | PARTIAL |
| Parallel PGGB | PASS | PASS | PASS | NOT_VERIFIED | NOT_RUN | PARTIAL |
| DNAnexus applets | PASS | PASS | NOT_RUN | NOT_VERIFIED | NOT_RUN | IMPLEMENTED_NOT_VERIFIED |
| Parallel timing | PASS | PASS | NOT_RUN | NOT_VERIFIED | NOT_RUN | IMPLEMENTED_NOT_VERIFIED |
| Stitch | PASS | PASS | PASS | NOT_RUN | NOT_RUN | **SYNTHETIC_VALIDATED** |
| Boundary validation | PASS | PASS | PASS (2/2) | NOT_RUN | NOT_RUN | SYNTHETIC_VALIDATED |
| Path equivalence | PASS | PASS | PASS (EQUIVALENT) | NOT_RUN | NOT_RUN | SYNTHETIC_VALIDATED |
| Variant equivalence | PASS | PASS | NOT_RUN | NOT_RUN | NOT_RUN | NOT_RUN |
| Web dataset (A6) | PASS | PASS | PASS | N/A | N/A | PASS |
| Graph explorer | PASS | PASS | PASS | N/A | N/A | PASS |
| Dashboard | PASS | PASS | PASS | N/A | N/A | PASS |
| Chunks | PASS | PASS | PASS | N/A | N/A | PASS |
| Compare | PARTIAL (stub) | PASS | PASS | N/A | N/A | PARTIAL |
| Vercel readiness | PASS | PASS | PASS | N/A | N/A | PASS |

## What Was Fixed In This Audit

1. **`scripts/run_pggb.py`** — Removed stale/invalid `-w` flag, wrong param names
   (`kmer_length`, `window_size`, `map_pct_id`, `noise_filter`), and `:latest`
   fallback. Now imports `gen_pggb_config.load_config()` and uses canonical
   `-p -s -K -k -j -e` flags + pinned digest + `samtools faidx`.

2. **`dnanexus/docker_helper.sh`** — `run_pggb()` removed undefined variables
   (`$MIN_ID`, `$KMER`, `$WINDOW`, `$MAP_PCT`, `$NOISE`) and the `-w` flag.
   Now reads canonical config JSON.

3. **`dnanexus/run_pipeline.sh`** — Updated stale "stitch NOT_IMPLEMENTED"
   messages to reflect `overlap_aware` stitch.

4. **`Makefile`** — Removed `|| true` masking scientific benchmark failure.

5. **`scripts/sync_web_results.py`** — Added `.bam`, `.cram`, `.fasta.gz` to
   forbidden genomic extensions.

6. **`pipeline/benchmark/benchmark_variants.sh`** — Removed `|| true` masking
   vg/truvari failures; added pinned Docker vg image fallback.

7. **`pipeline/parallel/build_chunk.sh`** — Removed hardcoded PGGB flags and
   `:latest` fallback. Now delegates to `scripts/run_pggb.py`.

## Stitch Status

**SYNTHETIC_VALIDATED** — Quang's overlap-aware stitch is merged into main
(commit `d64d3bb`, PR #1 + #3). `make demo` produces:
- Baseline: 502n 497e
- Merged: 511n 506e, 5 paths, 5 components
- Verdict: EQUIVALENT
- Boundaries: 2/2 PASS

This is **NOT** real-HPRC validated. Equivalence is a synthetic structural
comparison only.

## DNAnexus Status

**NOT_VERIFIED** — DNAnexus auth token expired (inactivity timeout) at audit
time. Applet IDs previously built were not re-verified.

## Remaining Blockers

1. **Real HPRC validation** — requires DNAnexus re-auth + full 1 Mb run.
2. **Variant equivalence** — requires real baseline + real stitched graph +
   vg/truvari on DNAnexus.

## DO NOT CLAIM

- Stitch is real-HPRC validated (only synthetic).
- Any real equivalence/benchmark numbers (none exist yet).
- Graph branches are SVs (unvalidated).
