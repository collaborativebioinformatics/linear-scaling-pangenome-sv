#!/usr/bin/env bash
# download_results.sh — Download benchmark results and web data from DNAnexus.
# Run LOCALLY (or inside workstation) after pipeline completes.
set -euo pipefail

echo "=== Download Results from DNAnexus ==="
echo ""

PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
if [ -z "$PROJECT_ID" ]; then
    echo "ERROR: DX_PROJECT_CONTEXT_ID unset. Run inside DNAnexus or set manually."
    exit 1
fi
echo "Project: $PROJECT_ID"
echo ""

# Benchmark report
echo "Downloading benchmark report..."
dx download "$PROJECT_ID:/benchmark/report.json" -o results/benchmark/report.json 2>/dev/null \
    && echo "  -> results/benchmark/report.json" \
    || echo "  (not yet available)"

# Merged graph
echo "Downloading merged graph..."
dx download "$PROJECT_ID:/graphs/merged/merged.gfa" -o results/merge/merged.gfa 2>/dev/null \
    && echo "  -> results/merge/merged.gfa" \
    || echo "  (not yet available)"

# Web data
echo "Downloading web data..."
dx download "$PROJECT_ID:/web/latest.json" -o web/public/data/latest.json 2>/dev/null \
    && echo "  -> web/public/data/latest.json" \
    || echo "  (not yet available)"

echo ""
echo "=== Download Complete ==="