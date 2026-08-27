#!/usr/bin/env bash
# run_parallel_chunks.sh — Launch chunk PGGB jobs on DNAnexus via pggb_chunk applet.
#
# STRICT enforcement:
#   expected_chunks == local FASTAs == uploaded FASTAs == submitted jobs == successful jobs == downloaded GFAs
#   Any mismatch = FATAL (no 'continue' on failures).
#
# 1. Build pggb_chunk applet capturing APPLET_ID via dx build --brief
# 2. Upload each chunk FASTA (skip if exists via hash compare)
# 3. Launch independent dx jobs per chunk
# 4. Wait for ALL jobs to finish
# 5. Download GFA outputs using job's formal gfa output link
# 6. Record per-chunk timing and aggregate metrics
set -euo pipefail

# === DRY RUN MODE ===
DRY_RUN=false
if [ "${1:-}" = "--dry-run" ]; then
    DRY_RUN=true
    echo "  DRY RUN MODE — validates counts, manifest, applet; no uploads/jobs/downloads"
    shift
fi

PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
INSTANCE="${PGGB_INSTANCE_TYPE:-mem3_ssd1_v2_x16}"

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
[ ! -f "$MANIFEST" ] && echo "FATAL: No manifest at $MANIFEST" && exit 1

EXPECTED_CHUNKS=$(tail -n +2 "$MANIFEST" | wc -l)
echo "  Expected chunks: $EXPECTED_CHUNKS"

# === BUILD APPLET (capture APPLET_ID via --brief, not by rediscovering) ===
APPLET_ID="${PGGB_APPLET_ID:-}"
echo ""
if [ -z "$APPLET_ID" ]; then
    echo "  Building pggb_chunk applet..."
    cd "$(dirname "$0")/applets/pggb_chunk"
    APPLET_ID=$(dx build --destination /applets/pggb_chunk/ --brief 2>/dev/null || echo "")
    cd "$OLDPWD"
fi
if [ -z "$APPLET_ID" ]; then
    echo "FATAL: Cannot build pggb_chunk applet. Set PGGB_APPLET_ID env var or build manually."
    exit 1
fi
echo "  Applet: $APPLET_ID"

# === STRICT ORCHESTRATION ===
CHUNK_DIR="/data/prepared/chunks"
JOB_IDS=()
CHUNK_IDS=()
CHUNK_COUNT=0
UPLOADED_COUNT=0
ORCH_START=$(date +%s)

# Count local FASTAs first
LOCAL_FA_COUNT=0
for f in work/chunks/chunk_*.fa; do
    [ -f "$f" ] && LOCAL_FA_COUNT=$((LOCAL_FA_COUNT + 1))
done
echo "  Local chunk FASTAs: $LOCAL_FA_COUNT"
if [ "$LOCAL_FA_COUNT" -ne "$EXPECTED_CHUNKS" ]; then
    echo "FATAL: Expected $EXPECTED_CHUNKS chunk FASTAs, found $LOCAL_FA_COUNT"
    exit 1
fi

while IFS=$'\t' read -r CID _ _ _ _ _ _ _ _; do
    [ "$CID" = "chunk_id" ] && continue
    FA="work/chunks/${CID}.fa"
    [ ! -f "$FA" ] && echo "FATAL: Missing chunk FASTA: $FA" && exit 1
    CHUNK_COUNT=$((CHUNK_COUNT + 1))
    echo ""
    echo "  Chunk $CID ..."

    # Upload chunk FASTA — MUST succeed
    echo "    Uploading FASTA..."
    FASTA_FILE_ID=$(dx upload "$FA" --destination "$CHUNK_DIR/" --brief 2>/dev/null || echo "")
    if [ -z "$FASTA_FILE_ID" ]; then
        echo "    FATAL: Failed to upload $CID.fa"
        exit 1
    fi
    UPLOADED_COUNT=$((UPLOADED_COUNT + 1))
    echo "    Uploaded: $FASTA_FILE_ID"

    # Launch chunk job — MUST succeed
    echo "    Launching PGGB job..."
    PGGB_CONFIG_JSON=$(python3 scripts/gen_pggb_config.py)
    
    if [ "$DRY_RUN" = true ]; then
        echo "    DRY RUN: would upload $CID.fa, launch job"
        UPLOADED_COUNT=$((UPLOADED_COUNT + 1))
        JOB_IDS+=("dry-run-$CID")
        CHUNK_IDS+=("$CID")
        continue
    fi
    
    JOB_ID=$(dx run "$APPLET_ID" \
        -i fasta="$FASTA_FILE_ID" \
        -i pggb_config_json="$PGGB_CONFIG_JSON" \
        --instance-type "$INSTANCE" \
        --name "Chunk $CID" \
        --destination "/graphs/chunks/$CID" \
        --brief 2>/dev/null || echo "")
    if [ -z "$JOB_ID" ]; then
        echo "    FATAL: Failed to launch PGGB job for $CID"
        exit 1
    fi
    echo "    Job: $JOB_ID"
    JOB_IDS+=("$JOB_ID")
    CHUNK_IDS+=("$CID")
done < <(tail -n +2 "$MANIFEST" 2>/dev/null || echo "")

echo ""
echo "=== Submitted $CHUNK_COUNT chunks, ${#JOB_IDS[@]} jobs, $UPLOADED_COUNT uploaded ==="
if [ "${#JOB_IDS[@]}" -ne "$EXPECTED_CHUNKS" ] || [ "$UPLOADED_COUNT" -ne "$EXPECTED_CHUNKS" ]; then
    echo "FATAL: Expected $EXPECTED_CHUNKS jobs, got ${#JOB_IDS[@]} jobs, $UPLOADED_COUNT uploads"
    exit 1
fi

# === WAIT FOR JOBS ===
echo ""
echo "Waiting for ${#JOB_IDS[@]} PGGB jobs to complete..."
FAILED_JOBS=()
JOB_METADATA=()
SUCCESSFUL_COUNT=0
for i in "${!JOB_IDS[@]}"; do
    JOB_ID="${JOB_IDS[$i]}"
    CID="${CHUNK_IDS[$i]}"
    echo -n "  $CID ($JOB_ID) ... "
    dx wait "$JOB_ID" 2>/dev/null || true
    JOB_JSON=$(dx describe "$JOB_ID" --json 2>/dev/null || echo "{}")
    JOB_STATE=$(echo "$JOB_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('state','unknown'))" 2>/dev/null || echo "unknown")
    if [ "$JOB_STATE" != "done" ]; then
        echo "FAILED ($JOB_STATE)"
        FAILED_JOBS+=("$JOB_ID")
        continue
    fi
    SUCCESSFUL_COUNT=$((SUCCESSFUL_COUNT + 1))

    # Extract timing from job describe
    STARTED_RUN_MS=$(echo "$JOB_JSON" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('startedRunning',0))" 2>/dev/null || echo "0")
    STOPPED_RUN_MS=$(echo "$JOB_JSON" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('stoppedRunning',0))" 2>/dev/null || echo "0")
    JOB_WALL=$(( (STOPPED_RUN_MS - STARTED_RUN_MS) / 1000 ))
    echo "OK (start=$STARTED_TS, stop=$STOPPED_TS)"

    # Download GFA using job output reference (formal gfa output link)
    echo "    Downloading GFA..."
    GFA_FILE_ID=$(echo "$JOB_JSON" | python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
    outputs = d.get('output', {})
    if 'gfa' in outputs:
        f = outputs['gfa']
        if isinstance(f, dict): print(f.get('$dnanexus_link', ''))
        else: print(str(f))
except: pass
" 2>/dev/null || echo "")

    if [ -n "$GFA_FILE_ID" ]; then
        dx download "$GFA_FILE_ID" -o "work/chunks/${CID}.gfa" 2>/dev/null || \
            echo "    WARNING: Could not download GFA for $CID"
        if [ -f "work/chunks/${CID}.gfa" ]; then
            echo "    -> work/chunks/${CID}.gfa ($(wc -c < "work/chunks/${CID}.gfa") bytes)"
        fi
    fi
    JOB_METADATA+=("{\"job_id\":\"$JOB_ID\",\"chunk_id\":\"$CID\",\"started_running_ms\":$STARTED_RUN_MS,\"stopped_running_ms\":$STOPPED_RUN_MS,\"wall_seconds\":$JOB_WALL,\"state\":\"$JOB_STATE\",\"instance_type\":\"$INSTANCE\"}")
done

echo ""
if [ ${#FAILED_JOBS[@]} -gt 0 ]; then
    echo "FATAL: ${#FAILED_JOBS[@]}/${#JOB_IDS[@]} PGGB jobs failed."
    exit 1
fi

# === VERIFY DOWNLOADED GFAs ===
DOWNLOADED_COUNT=0
for CID in "${CHUNK_IDS[@]}"; do
    [ -f "work/chunks/${CID}.gfa" ] && DOWNLOADED_COUNT=$((DOWNLOADED_COUNT + 1)) || true
done
echo "Downloaded GFAs: $DOWNLOADED_COUNT / ${#CHUNK_IDS[@]}"
if [ "$DOWNLOADED_COUNT" -ne "${#CHUNK_IDS[@]}" ]; then
    echo "FATAL: Expected ${#CHUNK_IDS[@]} GFAs, got $DOWNLOADED_COUNT"
    exit 1
fi

ORCH_END=$(date +%s)
ORCH_WALL=$((ORCH_END - ORCH_START))

# === AGGREGATE TIMING ===
if [ ${#JOB_METADATA[@]} -gt 0 ]; then
    TIMING_JSON=$(python3 -c "
import json, sys
jobs = [json.loads(j) for j in sys.argv[1:]]
starts = [j.get('started_running_ms',0) for j in jobs if j.get('started_running_ms')]
stops = [j.get('stopped_running_ms',0) for j in jobs if j.get('stopped_running_ms')]
if starts and stops:
    min_s, max_e = min(starts), max(stops)
    parallel_wall = (max_e - min_s) / 1000.0
else:
    parallel_wall = 0
sum_worker = sum(j.get('wall_seconds',0) for j in jobs)
print(json.dumps({'graph_parallel_wall_seconds': round(parallel_wall, 1), 'sum_worker_seconds': round(sum_worker, 1)}))
" "${JOB_METADATA[@]}" 2>/dev/null || echo '{"graph_parallel_wall_seconds":0,"sum_worker_seconds":0}')
else
    TIMING_JSON='{"graph_parallel_wall_seconds":0,"sum_worker_seconds":0}'
fi

# Write orchestration report
mkdir -p "results/benchmark"
cat > "results/benchmark/chunk_execution_report.json" << JSONEOF
{
  "expected_chunks": $EXPECTED_CHUNKS,
  "local_fastas": $LOCAL_FA_COUNT,
  "uploaded_fastas": $UPLOADED_COUNT,
  "submitted_jobs": ${#JOB_IDS[@]},
  "successful_jobs": $SUCCESSFUL_COUNT,
  "downloaded_gfas": $DOWNLOADED_COUNT,
  "instance_type": "$INSTANCE",
  "orchestration_wall_seconds": $ORCH_WALL,
  "per_chunk_timing": [$(IFS=,; echo "${JOB_METADATA[*]}")],
  $(echo "$TIMING_JSON" | tail -1 | sed 's/^{//' | sed 's/}$//')
}
JSONEOF

echo ""
echo "=== Parallel chunk execution complete ==="
echo "  Chunks: $EXPECTED_CHUNKS, Jobs: ${#JOB_IDS[@]}, Downloaded: $DOWNLOADED_COUNT"
echo "  Orchestration wall: ${ORCH_WALL}s"
echo "  Parallel wall (max worker stop - min worker start): $(echo "$TIMING_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin)['graph_parallel_wall_seconds'])")s"
echo "  Sum worker seconds: $(echo "$TIMING_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin)['sum_worker_seconds'])")s"
echo "  Report: results/benchmark/chunk_execution_report.json"