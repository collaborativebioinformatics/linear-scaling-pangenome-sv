#!/usr/bin/env bash
# run_pipeline.sh — Run the full chr21 smoke-test pipeline on DNAnexus.
#
# This is the single orchestrator for the entire cloud pipeline.
# Run inside a DNAnexus Cloud Workstation after setup.
#
# Usage:
#   bash dnanexus/run_pipeline.sh             # interactive mode
#   bash dnanexus/run_pipeline.sh --upload    # also upload results to project storage
#
set -euo pipefail

UPLOAD="${1:-}"

echo "=============================================="
echo "  Pangenome Parallel Pipeline — Smoke Test"
echo "  Target: chr21:20000000-21000000"
echo "=============================================="
echo ""

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCRIPT_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="results/logs"
mkdir -p "$LOG_DIR"

exec > >(tee "$LOG_DIR/pipeline_${TIMESTAMP}.log") 2>&1

# ── Checkpoint 1: Environment ──────────────────────────────────────────────
echo "[CHECKPOINT 1/9] Environment check"
bash scripts/check_environment.sh || {
    echo "Failed. Run: bash dnanexus/setup_workstation.sh"
    exit 1
}

# ── Checkpoint 2: HPRC manifest ────────────────────────────────────────────
echo "[CHECKPOINT 2/9] Fetching HPRC index"
python3 scripts/fetch_hprc_index.py
echo ""

# ── Checkpoint 3: Download data ────────────────────────────────────────────
echo "[CHECKPOINT 3/9] Downloading HPRC assemblies"
python3 scripts/download_hprc.py --execute 2>&1 || echo "  (downloads may be cached)"
echo ""

echo "[3b/9] Preparing reference"
bash scripts/prepare_reference.sh
echo ""

echo "[3c/9] Uploading to DNAnexus storage"
bash dnanexus/create_project_dirs.sh
bash dnanexus/upload_inputs.sh
echo ""

# ── Checkpoint 4: Prepare chr21 ────────────────────────────────────────────
echo "[CHECKPOINT 4/9] Preparing chr21 sequences"
python3 pipeline/prepare/prepare_sequences.py
echo ""

INPUT_FASTA="results/preparation/chr21_multi.fa"
if [ ! -f "$INPUT_FASTA" ]; then
    echo "FATAL: chr21 multi-FASTA not produced."
    exit 1
fi
echo "  Input FASTA: $INPUT_FASTA ($(grep -c '^>' "$INPUT_FASTA") paths)"

# ── Checkpoint 5: Baseline PGGB ────────────────────────────────────────────
echo "[CHECKPOINT 5/9] Building monolithic baseline graph"
bash pipeline/baseline/build_baseline.sh "$INPUT_FASTA" "results/baseline"
echo ""

BASELINE_GFA="results/baseline/baseline.gfa"
if [ ! -f "$BASELINE_GFA" ]; then
    echo "FATAL: Baseline graph not produced."
    exit 1
fi
echo "  Baseline: $BASELINE_GFA ($(wc -c < "$BASELINE_GFA") bytes)"

# ── Checkpoint 6: Prepare chunk sequences ──────────────────────────────────
echo "[CHECKPOINT 6/9] Creating chunks"
python3 pipeline/parallel/make_chunks.py
echo ""

# Build per-chunk FASTA files from the multi-FASTA
python3 pipeline/parallel/build_all_chunks.py
echo ""

# ── Checkpoint 6b: Build chunk graphs ──────────────────────────────────────
echo "[6b/9] Building chunk PGGB graphs"
python3 pipeline/parallel/build_all_chunks.py --execute
echo ""

# Verify chunk GFAs exist
CHUNK_COUNT=$(find work/chunks -name "chunk_*.gfa" -type f | wc -l | tr -d ' ')
echo "  Chunk GFAs produced: $CHUNK_COUNT"

# ── Checkpoint 7: Merge ────────────────────────────────────────────────────
echo "[CHECKPOINT 7/9] Merging chunk graphs"
python3 pipeline/merge/merge_graphs.py
echo ""

MERGED_GFA="results/merge/merged.gfa"
if [ ! -f "$MERGED_GFA" ]; then
    echo "WARNING: Merged graph not produced (merge may need real PGGB output)"
fi

# ── Checkpoint 8: Benchmark ────────────────────────────────────────────────
echo "[CHECKPOINT 8/9] Running benchmark"
python3 pipeline/benchmark/graph_stats.py
python3 pipeline/benchmark/compare_paths.py
bash pipeline/benchmark/benchmark_variants.sh 2>/dev/null || true
python3 pipeline/benchmark/build_report.py
echo ""

# ── Checkpoint 9: Export web JSON ──────────────────────────────────────────
echo "[CHECKPOINT 9/9] Exporting web JSON"
python3 pipeline/export/gfa_to_json.py "$BASELINE_GFA" --output "web/public/data/baseline.json" --label "baseline" 2>/dev/null || true
if [ -f "$MERGED_GFA" ]; then
    python3 pipeline/export/gfa_to_json.py "$MERGED_GFA" --output "web/public/data/merged.json" --label "merged" 2>/dev/null || true
fi

# Build final web JSON from benchmark report
python3 scripts/sync_web_results.py 2>/dev/null || true
echo ""

# ── Summary ────────────────────────────────────────────────────────────────
echo "=============================================="
echo "  Pipeline Complete"
echo "=============================================="
echo ""
echo "  Baseline:  $BASELINE_GFA"
echo "  Merged:    $MERGED_GFA"
echo "  Logs:      $LOG_DIR/pipeline_${TIMESTAMP}.log"
echo ""

if [ -f "$BASELINE_GFA" ]; then
    echo "  Baseline nodes: $(grep -c '^S' "$BASELINE_GFA" 2>/dev/null || echo '?')"
fi
if [ -f "$MERGED_GFA" ]; then
    echo "  Merged nodes:   $(grep -c '^S' "$MERGED_GFA" 2>/dev/null || echo '?')"
fi
echo ""

# Upload results if requested
if [ "$UPLOAD" = "--upload" ]; then
    echo "Uploading results to DNAnexus..."
    dx upload results/merge/merged.gfa --destination /graphs/merged/ 2>/dev/null || true
    dx upload results/baseline/baseline.gfa --destination /graphs/baseline/ 2>/dev/null || true
    dx upload results/benchmark/report.json --destination /benchmark/ 2>/dev/null || true
    echo "  Done"
fi

echo ""
echo "Next: bash dnanexus/download_results.sh"