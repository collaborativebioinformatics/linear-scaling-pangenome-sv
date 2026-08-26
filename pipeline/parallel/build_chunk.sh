#!/usr/bin/env bash
# build_chunk.sh — Build PGGB graph for a single chunk.
# REAL-DATA STEP: Requires pggb and chunk FASTA.
# TODO (Michael): Verify PGGB params work at chunk scale.
# TODO (Quang): Add container fallback.
set -euo pipefail

CHUNK_ID="${1:?Usage: build_chunk.sh <chunk_id>}"
INPUT="work/chunks/${CHUNK_ID}.fa"
OUTDIR="work/chunks/${CHUNK_ID}"
THREADS="${PGGB_THREADS:-16}"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Chunk FASTA not found: $INPUT"
    exit 1
fi

mkdir -p "$OUTDIR" "results/logs"
echo "Building chunk: $CHUNK_ID (threads=$THREADS)"

START=$(date +%s)

pggb \
    -i "$INPUT" \
    -o "$OUTDIR" \
    -t "$THREADS" \
    -n "$(grep -c '^>' "$INPUT")" \
    -p 90 \
    -s 5000 \
    -k 29 \
    -j 0 \
    -e 0 \
    2>&1 | tee "results/logs/chunk_${CHUNK_ID}.log"

END=$(date +%s)
echo "$CHUNK_ID finished in $((END - START))s"

GFA_FILE=$(find "$OUTDIR" -name "*.gfa" -type f | head -1)
if [ -n "$GFA_FILE" ]; then
    cp "$GFA_FILE" "work/chunks/${CHUNK_ID}.gfa"
    echo "  GFA: work/chunks/${CHUNK_ID}.gfa"
fi