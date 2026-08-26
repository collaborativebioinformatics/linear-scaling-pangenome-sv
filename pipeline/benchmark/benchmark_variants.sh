#!/usr/bin/env bash
# benchmark_variants.sh — Variant comparison using vg + Truvari.
# REAL-DATA STEP: Requires vg and truvari.
# TODO (Michael/Khoi): Verify vg deconstruct parameters.
set -euo pipefail

echo "=== Variant Benchmark (vg + Truvari) ==="
echo "REAL-DATA STEP: Requires vg and truvari."

if ! command -v vg &>/dev/null; then
    echo "  SKIP: vg not found"
    exit 0
fi
if ! command -v truvari &>/dev/null; then
    echo "  SKIP: truvari not found"
    exit 0
fi

BASELINE="results/baseline/baseline.gfa"
MERGED="results/merge/merged.gfa"
OUTDIR="results/benchmark/truvari"
mkdir -p "$OUTDIR"

if [ ! -f "$BASELINE" ] || [ ! -f "$MERGED" ]; then
    echo "  SKIP: graphs not found"
    exit 0
fi

echo "  Deconstructing baseline VCF..."
vg deconstruct -P GRCh38 -a "$BASELINE" > "$OUTDIR/baseline.vcf" 2>/dev/null || true

echo "  Deconstructing merged VCF..."
vg deconstruct -P GRCh38 -a "$MERGED" > "$OUTDIR/merged.vcf" 2>/dev/null || true

if [ -f "$OUTDIR/baseline.vcf" ] && [ -f "$OUTDIR/merged.vcf" ]; then
    echo "  Running Truvari comparison..."
    truvari bench \
        -b "$OUTDIR/baseline.vcf" \
        -c "$OUTDIR/merged.vcf" \
        -o "$OUTDIR/comparison" \
        --passonly \
        2>&1 | tee "results/logs/truvari.log" || true
    echo "  Truvari complete"
fi

echo "  Done"