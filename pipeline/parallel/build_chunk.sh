#!/usr/bin/env bash
# build_chunk.sh — Build PGGB graph for a single chunk (LOCAL / development).
# LEGACY/DEVELOPMENT_ONLY: the canonical chunk execution path is the DNAnexus
# pggb_chunk applet. This local wrapper delegates to scripts/run_pggb.py so it
# uses the SAME canonical config/pipeline.yaml parameters (no hard-coded flags).
# Usage: bash pipeline/parallel/build_chunk.sh <chunk_id>
set -euo pipefail

CHUNK_ID="${1:?Usage: build_chunk.sh <chunk_id>}"
INPUT="work/chunks/${CHUNK_ID}.fa"
OUTDIR="work/chunks/${CHUNK_ID}"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Chunk FASTA not found: $INPUT"
    exit 1
fi

mkdir -p "$OUTDIR" "results/logs"

echo "Building chunk: $CHUNK_ID (canonical config from config/pipeline.yaml)"

START=$(date +%s)
python3 scripts/run_pggb.py "$INPUT" "$OUTDIR"
END=$(date +%s)
echo "$CHUNK_ID finished in $((END - START))s"

if [ -f "$OUTDIR/final.gfa" ]; then
    cp "$OUTDIR/final.gfa" "work/chunks/${CHUNK_ID}.gfa"
    echo "  GFA: work/chunks/${CHUNK_ID}.gfa"
else
    echo "  WARNING: No final.gfa produced for $CHUNK_ID"
    exit 1
fi