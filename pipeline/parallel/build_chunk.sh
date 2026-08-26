#!/usr/bin/env bash
# build_chunk.sh — Build PGGB graph for a single chunk via Docker container.
# Usage: bash pipeline/parallel/build_chunk.sh <chunk_id>

set -euo pipefail

CHUNK_ID="${1:?Usage: build_chunk.sh <chunk_id>}"
INPUT="work/chunks/${CHUNK_ID}.fa"
OUTDIR="work/chunks/${CHUNK_ID}"
THREADS="${PGGB_THREADS:-16}"
PGGB_IMAGE="${PGGB_IMAGE:-ghcr.io/pangenome/pggb:latest}"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Chunk FASTA not found: $INPUT"
    exit 1
fi

mkdir -p "$OUTDIR" "results/logs"

if ! docker image inspect "$PGGB_IMAGE" &>/dev/null; then
    echo "Pulling PGGB container..."
    docker pull "$PGGB_IMAGE"
fi

INPUT_ABS="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
OUTDIR_ABS="$(cd "$OUTDIR" && pwd)"

echo "Building chunk: $CHUNK_ID (threads=$THREADS)"

START=$(date +%s)

docker run --rm \
    -v "$(dirname "$INPUT_ABS")":/data/input:ro \
    -v "$OUTDIR_ABS":/data/output \
    "$PGGB_IMAGE" \
    pggb \
        -i "/data/input/$(basename "$INPUT")" \
        -o "/data/output" \
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
else
    echo "  WARNING: No GFA produced for $CHUNK_ID"
fi