#!/usr/bin/env bash
# build_baseline.sh — Build monolithic PGGB graph for the smoke-test interval.
#
# Runs PGGB via Docker container. Requires:
#   - Multi-haplotype FASTA at INPUT path
#   - Docker with ghcr.io/pangenome/pggb:latest pulled
#
# Usage:
#   bash pipeline/baseline/build_baseline.sh [input.fa] [output_dir] [threads] [ref_name]

set -euo pipefail

INPUT="${1:-results/preparation/chr21_multi.fa}"
OUTDIR="${2:-results/baseline}"
THREADS="${3:-16}"
REF="${4:-GRCh38}"
PGGB_IMAGE="${PGGB_IMAGE:-ghcr.io/pangenome/pggb:latest}"

mkdir -p "$OUTDIR" "results/logs"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input FASTA not found: $INPUT"
    echo "Run preparation step first."
    exit 1
fi

# PGGB/wfmash require a bgzip-compressed, faidx-indexed FASTA.
if [[ "$INPUT" != *.gz ]]; then
    if [ -f "${INPUT}.gz" ]; then
        INPUT="${INPUT}.gz"
    else
        echo "Compressing input with bgzip..."
        bgzip -@ "$THREADS" -k "$INPUT"
        INPUT="${INPUT}.gz"
    fi
fi
if [ ! -f "${INPUT}.fai" ] || [ ! -f "${INPUT}.gzi" ]; then
    echo "Indexing $INPUT ..."
    samtools faidx "$INPUT"
fi
NSEQ=$(wc -l < "${INPUT}.fai" | tr -d ' ')
echo "Input sequences: $NSEQ"

if ! docker image inspect "$PGGB_IMAGE" &>/dev/null; then
    echo "Pulling PGGB container: $PGGB_IMAGE"
    docker pull "$PGGB_IMAGE"
fi

INPUT_ABS="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
OUTDIR_ABS="$(cd "$(dirname "$OUTDIR")" && pwd)/$(basename "$OUTDIR")"
RESULTS_ABS="$(cd "$(dirname "results/logs")/.." && pwd)/logs"

echo "=== Baseline PGGB Graph ==="
echo "Input: $INPUT"
echo "Output: $OUTDIR"
echo "Threads: $THREADS"
echo "Container: $PGGB_IMAGE"

START=$(date +%s)

docker run --rm \
    -v "$(dirname "$INPUT_ABS")":/data/input:ro \
    -v "$OUTDIR_ABS":/data/output \
    -v "$RESULTS_ABS":/data/logs \
    "$PGGB_IMAGE" \
    pggb \
        -i "/data/input/$(basename "$INPUT")" \
        -o "/data/output" \
        -t "$THREADS" \
        -n "$NSEQ" \
        -p 90 \
        -s 5000 \
        -k 29 \
        -w 50000 \
        -j 0 \
        -e 0 \
        2>&1 | tee "results/logs/baseline.log"

END=$(date +%s)
DURATION=$((END - START))
echo "PGGB finished in ${DURATION}s"

# Locate the GFA output (PGGB may produce subdirectories)
GFA_FILE=$(find "$OUTDIR" -name "*.gfa" -type f 2>/dev/null | head -1)
if [ -z "$GFA_FILE" ]; then
    echo "ERROR: No GFA produced by PGGB"
    echo "Check: ls -la $OUTDIR/"
    ls -la "$OUTDIR/" 2>/dev/null || true
    exit 1
fi

cp "$GFA_FILE" "$OUTDIR/baseline.gfa"
echo "Canonical: $OUTDIR/baseline.gfa"

cat > "$OUTDIR/run_metadata.json" << JSONEOF
{
  "method": "monolithic",
  "target": "chr21",
  "input_paths": $NSEQ,
  "threads": $THREADS,
  "wall_seconds": $DURATION,
  "peak_memory_kb": null,
  "status": "completed",
  "container": "$PGGB_IMAGE"
}
JSONEOF

echo "Baseline complete in ${DURATION}s"