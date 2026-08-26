#!/usr/bin/env bash
# run_smoke_test.sh — Run the full smoke test pipeline.
# TODO (Ali/Khoi): Add error handling and logging per stage.
# TODO (Quang): Parameterize chunk workflow.
set -euo pipefail

echo "=== Smoke Test Pipeline ==="
echo "Target: chr21:20000000-21000000"
echo ""

# Step 1: Check environment
echo "[1/8] Checking environment..."
bash scripts/check_environment.sh

# Step 2: Fetch HPRC index
echo "[2/8] Fetching HPRC index..."
python3 scripts/fetch_hprc_index.py

# Step 3: Download assemblies
echo "[3/8] Downloading assemblies..."
python3 scripts/download_hprc.py

# Step 4: Prepare chr21
echo "[4/8] Preparing chr21..."
python3 pipeline/prepare/prepare_sequences.py

# Step 5: Baseline PGGB graph
echo "[5/8] Building baseline graph..."
bash pipeline/baseline/build_baseline.sh

# Step 6: Parallel chunks
echo "[6/8] Creating chunks..."
python3 pipeline/parallel/make_chunks.py
python3 pipeline/parallel/build_all_chunks.py

# Step 7: Merge
echo "[7/8] Merging graphs..."
python3 pipeline/merge/merge_graphs.py

# Step 8: Benchmark
echo "[8/8] Running benchmark..."
python3 pipeline/benchmark/graph_stats.py
python3 pipeline/benchmark/compare_paths.py
python3 pipeline/benchmark/build_report.py

echo ""
echo "=== Smoke Test Complete ==="
echo "Results in results/"
echo "Run: python3 scripts/sync_web_results.py"