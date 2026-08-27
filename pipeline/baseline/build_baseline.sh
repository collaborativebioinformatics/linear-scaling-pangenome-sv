#!/usr/bin/env bash
# build_baseline.sh — Build monolithic PGGB graph for the smoke-test interval.
#
# Reads configuration from config/pipeline.yaml (single source of truth).
# Uses scripts/run_pggb.py for canonical PGGB execution.
# Outputs exactly one *final.gfa (FATAL if zero or >1).
#
# Usage:
#   bash pipeline/baseline/build_baseline.sh [input.fa] [output_dir]

set -euo pipefail

INPUT="${1:-results/preparation/chr21_multi.fa}"
OUTDIR="${2:-results/baseline}"
mkdir -p "$OUTDIR" "results/logs"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input FASTA not found: $INPUT"
    echo "Run preparation step first."
    exit 1
fi

echo "=== Baseline PGGB Graph ==="
echo "Input: $INPUT"
echo "Output: $OUTDIR"

START=$(date +%s)
python3 scripts/run_pggb.py "$INPUT" "$OUTDIR"
END=$(date +%s)
DURATION=$((END - START))

# Canonical output name
FINAL_GFA="$OUTDIR/final.gfa"
if [ -f "$FINAL_GFA" ]; then
    cp "$FINAL_GFA" "$OUTDIR/baseline.gfa"
    echo "Canonical: $OUTDIR/baseline.gfa"
else
    echo "FATAL: No final.gfa produced"
    exit 1
fi

echo "Baseline complete in ${DURATION}s"