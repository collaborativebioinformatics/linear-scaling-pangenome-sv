#!/usr/bin/env bash
# benchmark_variants.sh — Variant comparison using vg + Truvari.
# REAL-DATA STEP: Requires vg and truvari.
# When variant equivalence is explicitly enabled, failures are FATAL
# (no `|| true` masking of scientific benchmark failures).
set -euo pipefail

echo "=== Variant Benchmark (vg + Truvari) ==="
echo "REAL-DATA STEP: Requires vg and truvari."

# Use pinned Docker vg image when host vg is unavailable.
VG_BIN="vg"
if ! command -v vg &>/dev/null; then
    if command -v docker &>/dev/null; then
        VG_IMAGE="${VG_IMAGE:-quay.io/vgteam/vg:v1.74.1}"
        echo "  host vg not found; using Docker image $VG_IMAGE"
        VG_BIN="docker run --rm -v \$PWD:/data $VG_IMAGE vg"
    else
        echo "  SKIP: vg not found and no Docker available"
        exit 0
    fi
fi
if ! command -v truvari &>/dev/null; then
    echo "  SKIP: truvari not found"
    exit 0
fi

BASELINE="results/baseline/baseline.gfa"
MERGED="results/merge/merged.gfa"
OUTDIR="results/benchmark/truvari"
mkdir -p "$OUTDIR" results/logs

if [ ! -f "$BASELINE" ] || [ ! -f "$MERGED" ]; then
    echo "  SKIP: graphs not found (baseline + real merged required)"
    exit 0
fi

echo "  Deconstructing baseline VCF..."
$VG_BIN deconstruct -P GRCh38 -a "$BASELINE" > "$OUTDIR/baseline.vcf"

echo "  Deconstructing merged VCF..."
$VG_BIN deconstruct -P GRCh38 -a "$MERGED" > "$OUTDIR/merged.vcf"

echo "  Running Truvari comparison..."
truvari bench \
    -b "$OUTDIR/baseline.vcf" \
    -c "$OUTDIR/merged.vcf" \
    -o "$OUTDIR/comparison" \
    --passonly \
    2>&1 | tee "results/logs/truvari.log"

echo "  Truvari complete -> $OUTDIR/comparison"
echo "  Done"