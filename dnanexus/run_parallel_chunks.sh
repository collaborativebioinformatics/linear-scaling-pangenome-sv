#!/usr/bin/env bash
# run_parallel_chunks.sh — Launch chunk PGGB jobs via DNAnexus applet or local Docker.
# Uses the pggb_chunk applet. Builds it first if not found.
# Falls back to local Docker execution when applet is unavailable.
set -euo pipefail

PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
INSTANCE="${DX_INSTANCE_TYPE:-mem3_ssd1_v2_x16}"

echo "=== Run Parallel Chunks ==="

if [ -z "$PROJECT_ID" ]; then
    echo "ERROR: Not inside DNAnexus. Fallback to local Docker."
fi

# Ensure chunk manifest and FASTA files exist
python3 pipeline/parallel/make_chunks.py 2>/dev/null || true
python3 pipeline/parallel/build_all_chunks.py

MANIFEST="work/chunks/chunk_manifest.tsv"
[ ! -f "$MANIFEST" ] && echo "No manifest." && exit 1

# Check if pggb_chunk applet exists on DNAnexus
APPLET_ID=""
if [ -n "$PROJECT_ID" ]; then
    APPLET_ID=$(dx find data --path /applets/pggb_chunk --brief 2>/dev/null | head -1 || echo "")
    if [ -z "$APPLET_ID" ]; then
        echo "Building pggb_chunk applet in project..."
        cd dnanexus/applets/pggb_chunk
        dx build --destination /applets/pggb_chunk/ --brief 2>/dev/null || true
        cd /home/dnanexus/pangenome-parallel 2>/dev/null || cd ../..
        APPLET_ID=$(dx find data --path /applets/pggb_chunk --brief 2>/dev/null | head -1 || echo "")
    fi
fi

echo "Instance: $INSTANCE"
echo "Applet: ${APPLET_ID:-none (using local Docker)}"

# Process each chunk
JOB_IDS=()
while IFS=$'\t' read -r CID _ _ _ _ _ _ _ _; do
    [ "$CID" = "chunk_id" ] && continue
    FA="work/chunks/${CID}.fa"
    [ ! -f "$FA" ] && echo "  MISSING $FA" && continue

    if [ -n "$APPLET_ID" ]; then
        echo "  Launching $CID via DNAnexus applet..."
        JOB_ID=$(dx run "$APPLET_ID" \
            -i fasta="$FA" \
            --instance-type "$INSTANCE" \
            --name "Chunk $CID" \
            --destination "/graphs/chunks/$CID" \
            --brief 2>/dev/null || echo "")
        if [ -n "$JOB_ID" ]; then
            echo "    Job: $JOB_ID"
            JOB_IDS+=("$JOB_ID")
        else
            echo "    FAILED"
        fi
    else
        echo "  Building $CID via local Docker..."
        bash pipeline/parallel/build_chunk.sh "$CID"
    fi
done < <(tail -n +2 "$MANIFEST" 2>/dev/null || echo "")

echo ""
echo "${#JOB_IDS[@]} jobs submitted."
if [ -n "$APPLET_ID" ]; then
    echo "Monitor: dx find jobs --name 'Chunk *'"
fi
echo "After completion, merge: python3 pipeline/merge/merge_graphs.py"