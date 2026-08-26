#!/usr/bin/env bash
# run_parallel_chunks.sh — Submit PGGB chunk jobs as independent DNAnexus runs.
# Scatter-gather: launches each chunk as a separate dx job, then merges.
set -euo pipefail

echo "=== Run Parallel Chunks on DNAnexus ==="
echo ""

PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
INSTANCE="${DX_INSTANCE_TYPE:-mem3_ssd1_v2_x16}"

if [ -z "$PROJECT_ID" ]; then
    echo "ERROR: Not running inside DNAnexus."
    exit 1
fi

MANIFEST="work/chunks/chunk_manifest.tsv"
if [ ! -f "$MANIFEST" ]; then
    echo "Creating chunk manifest..."
    python3 pipeline/parallel/make_chunks.py
fi

echo "Instance type: $INSTANCE"
echo "Project: $PROJECT_ID"
echo ""

# Read chunks from manifest
JOB_IDS=()
CHUNK_COUNT=0
while IFS=$'\t' read -r CHUNK_ID _ START END _ _ _ _ _; do
    if [ "$CHUNK_ID" = "chunk_id" ]; then continue; fi  # skip header
    CHUNK_COUNT=$((CHUNK_COUNT + 1))

    INPUT_FASTA="/data/prepared/${CHUNK_ID}.fa"
    OUTPUT_DIR="/graphs/chunks/${CHUNK_ID}"

    echo "Submitting $CHUNK_ID ($START-$END)..."
    JOB_ID=$(dx run pggb \
        -i "fasta=$INPUT_FASTA" \
        --instance-type "$INSTANCE" \
        --name "Chunk $CHUNK_ID" \
        --destination "$OUTPUT_DIR" \
        --brief 2>/dev/null || echo "")

    if [ -n "$JOB_ID" ]; then
        echo "  Job: $JOB_ID"
        JOB_IDS+=("$JOB_ID")
    else
        echo "  FAILED to submit $CHUNK_ID"
    fi
done < <(tail -n +2 "$MANIFEST" 2>/dev/null || echo "")

echo ""
echo "Submitted $CHUNK_COUNT chunk jobs."
echo "Monitor with: dx find jobs --name 'Chunk *'"
echo ""

# Merge step (after all chunks complete)
echo "After all jobs complete, merge with:"
echo "  python3 pipeline/merge/merge_graphs.py"
echo "  dx upload results/merge/merged.gfa --destination /graphs/merged/"
echo ""
echo "=== Submission Complete ==="