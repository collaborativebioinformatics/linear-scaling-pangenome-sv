#!/usr/bin/env bash
# run_parallel_chunks.sh — Launch chunk PGGB jobs on DNAnexus via pggb_chunk applet.
#
# 1. Build/find the pggb_chunk applet (FATAL if unavailable — no local fallback)
# 2. Upload each chunk FASTA to /data/prepared/chunks/ (skip if exists)
# 3. Launch independent dx jobs per chunk
# 4. Wait for ALL jobs to finish
# 5. Download GFA outputs to work/chunks/<chunk_id>.gfa
# 6. Record per-chunk timing
set -euo pipefail

PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
INSTANCE="${DX_INSTANCE_TYPE:-mem3_ssd1_v2_x16}"

echo "=== Run Parallel Chunks on DNAnexus ==="
if [ -z "$PROJECT_ID" ]; then
    echo "FATAL: Not inside DNAnexus (DX_PROJECT_CONTEXT_ID unset)."
    exit 1
fi
echo "  Project: $PROJECT_ID"
echo "  Instance: $INSTANCE"

# Ensure chunk FASTAs exist
python3 pipeline/parallel/make_chunks.py 2>/dev/null || true
python3 pipeline/parallel/build_all_chunks.py

MANIFEST="work/chunks/chunk_manifest.tsv"
[ ! -f "$MANIFEST" ] && echo "No manifest." && exit 1

# Build/find pggb_chunk applet (FATAL if unavailable)
APPLET_ID=""
echo ""
echo "  Locating pggb_chunk applet..."
APPLET_ID=$(dx find data --path /applets/pggb_chunk --brief 2>/dev/null | head -1 || echo "")
if [ -z "$APPLET_ID" ]; then
    echo "  Building pggb_chunk applet from source..."
    cd dnanexus/applets/pggb_chunk
    dx build --destination /applets/pggb_chunk/ --brief
    cd "$OLDPWD"
    APPLET_ID=$(dx find data --path /applets/pggb_chunk --brief 2>/dev/null | head -1 || echo "")
fi
if [ -z "$APPLET_ID" ]; then
    echo "FATAL: pggb_chunk applet could not be built/found. Cannot run parallel PGGB without it."
    exit 1
fi
echo "  Applet: $APPLET_ID"

# Process each chunk
CHUNK_DIR="/data/prepared/chunks"
JOB_IDS=()
CHUNK_COUNT=0
PARALLEL_START=$(date +%s)

while IFS=$'\t' read -r CID _ _ _ _ _ _ _ _; do
    [ "$CID" = "chunk_id" ] && continue
    FA="work/chunks/${CID}.fa"
    [ ! -f "$FA" ] && echo "  MISSING $FA" && continue
    CHUNK_COUNT=$((CHUNK_COUNT + 1))

    # Upload chunk FASTA if not already in project storage
    echo ""
    echo "  Chunk $CID..."

    EXISTING_FILE_ID=$(dx find data --name "${CID}.fa" --path "$PROJECT_ID:$CHUNK_DIR" --brief 2>/dev/null | head -1 || echo "")
    if [ -n "$EXISTING_FILE_ID" ]; then
        echo "    FASTA already in project: $EXISTING_FILE_ID"
        FASTA_FILE_ID="$EXISTING_FILE_ID"
    else
        echo "    Uploading FASTA..."
        FASTA_FILE_ID=$(dx upload "$FA" --destination "$CHUNK_DIR/" --brief 2>/dev/null || echo "")
        if [ -z "$FASTA_FILE_ID" ]; then
            echo "    FAILED to upload $CID.fa"
            continue
        fi
        echo "    Uploaded: $FASTA_FILE_ID"
    fi

    # Launch chunk job
    echo "    Launching PGGB job..."
    JOB_ID=$(dx run "$APPLET_ID" \
        -i fasta="$FASTA_FILE_ID" \
        --instance-type "$INSTANCE" \
        --name "Chunk $CID" \
        --destination "/graphs/chunks/$CID" \
        --brief 2>/dev/null || echo "")

    if [ -n "$JOB_ID" ]; then
        echo "    Job: $JOB_ID"
        JOB_IDS+=("$JOB_ID")
    else
        echo "    FAILED to launch $CID"
    fi
done < <(tail -n +2 "$MANIFEST" 2>/dev/null || echo "")

echo ""
echo "=== Submitted $CHUNK_COUNT chunks, ${#JOB_IDS[@]} jobs ==="
echo ""

if [ ${#JOB_IDS[@]} -eq 0 ]; then
    echo "FATAL: No jobs were submitted."
    exit 1
fi

# Wait for all jobs to complete
echo "Waiting for ${#JOB_IDS[@]} PGGB jobs to complete..."
FAILED_JOBS=()
for JOB_ID in "${JOB_IDS[@]}"; do
    echo -n "  $JOB_ID ... "
    JOB_INFO=$(dx wait "$JOB_ID" 2>/dev/null || echo "failed")
    JOB_STATUS=$(dx describe "$JOB_ID" --json 2>/dev/null | python3 -c \
        "import sys,json; print(json.load(sys.stdin).get('state','unknown'))" 2>/dev/null || echo "unknown")
    if [ "$JOB_STATUS" = "done" ]; then
        echo "OK ($JOB_STATUS)"
    else
        echo "FAILED ($JOB_STATUS)"
        FAILED_JOBS+=("$JOB_ID")
    fi
done

echo ""
if [ ${#FAILED_JOBS[@]} -gt 0 ]; then
    echo "FATAL: ${#FAILED_JOBS[@]}/${#JOB_IDS[@]} PGGB jobs failed."
    for JID in "${FAILED_JOBS[@]}"; do
        echo "  Failed: $JID"
    done
    exit 1
fi

PARALLEL_END=$(date +%s)
PARALLEL_WALL=$((PARALLEL_END - PARALLEL_START))
echo "All ${#JOB_IDS[@]} PGGB jobs completed successfully."
echo "Parallel wall time: ${PARALLEL_WALL}s"

# Download GFA outputs
echo ""
echo "Downloading GFA outputs..."
for JOB_ID in "${JOB_IDS[@]}"; do
    JOB_NAME=$(dx describe "$JOB_ID" --json 2>/dev/null | python3 -c \
        "import sys,json; print(json.load(sys.stdin).get('name','unknown'))" 2>/dev/null || echo "unknown")
    CID="${JOB_NAME#Chunk }"
    OUTPUT_DIR="/graphs/chunks/$CID"
    echo "  $CID ..."
    dx download "$PROJECT_ID:$OUTPUT_DIR/merged.gfa" -o "work/chunks/${CID}.gfa" 2>/dev/null || \
        dx download "$PROJECT_ID:$OUTPUT_DIR/*.gfa" -o "work/chunks/${CID}.gfa" 2>/dev/null || \
        echo "    WARNING: GFA not found for $CID"
    if [ -f "work/chunks/${CID}.gfa" ]; then
        echo "    -> work/chunks/${CID}.gfa ($(wc -c < "work/chunks/${CID}.gfa") bytes)"
    fi
done

echo ""
echo "=== Parallel chunk execution complete ==="
echo "  Chunks: $CHUNK_COUNT, Jobs: ${#JOB_IDS[@]}, Wall: ${PARALLEL_WALL}s"
echo "  Merge: python3 pipeline/merge/merge_graphs.py"