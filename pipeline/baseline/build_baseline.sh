#!/usr/bin/env bash
# build_baseline.sh — Build monolithic PGGB graph.
# REAL-DATA STEP: Requires pggb and prepared multi-haplotype FASTA.
# TODO (Michael): Validate PGGB parameters for HPRC chr21 dataset.
# TODO (Khoi): Add container-based fallback when pggb not installed locally.
set -euo pipefail

INPUT="${1:-results/preparation/chr21_multi.fa}"
OUTDIR="${2:-results/baseline}"
THREADS="${3:-16}"
REF="${4:-GRCh38}"

mkdir -p "$OUTDIR" "results/logs"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input FASTA not found: $INPUT"
    echo "Run preparation step first."
    exit 1
fi

echo "=== Baseline PGGB Graph ==="
echo "Input: $INPUT"
echo "Output: $OUTDIR"
echo "Threads: $THREADS"

START=$(date +%s)

pggb \
    -i "$INPUT" \
    -o "$OUTDIR" \
    -t "$THREADS" \
    -n "$(grep -c '^>' "$INPUT")" \
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

GFA_FILE=$(find "$OUTDIR" -name "*.gfa" -type f 2>/dev/null | head -1)
if [ -z "$GFA_FILE" ]; then
    echo "ERROR: No GFA produced by PGGB"
    exit 1
fi
cp "$GFA_FILE" "$OUTDIR/baseline.gfa"

cat > "$OUTDIR/run_metadata.json" << JSONEOF
{
  "method": "monolithic",
  "target": "chr21",
  "input_paths": $(grep -c '^>' "$INPUT"),
  "threads": $THREADS,
  "wall_seconds": $DURATION,
  "peak_memory_kb": null,
  "status": "completed"
}
JSONEOF

echo "Baseline graph: $OUTDIR/baseline.gfa"