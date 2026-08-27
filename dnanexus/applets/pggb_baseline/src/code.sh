#!/usr/bin/env bash
# pggb_baseline applet — runs PGGB via Docker inside a DNAnexus job.
# Shares EXACT PGGB parameters with pggb_chunk via pggb_config_json string input.
# Same instance type, same threads, same image.
set -e -o pipefail

main() {
    echo "=== PGGB Baseline Graph Builder (DNAnexus) ==="
    echo "Input FASTA: $fasta_name"
    echo ""

    dx download "$fasta" -o input.fa
    NUM_PATHS=$(grep -c '^>' input.fa)
    echo "Paths in input: $NUM_PATHS"

    # Parse PGGB config from JSON string input (stdlib json only)
    CFG_JSON="$pggb_config_json"
    THREADS=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('threads',8))")
    MIN_ID=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('minimum_identity',90))")
    SEG_LEN=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('segment_length',5000))")
    MATCH_LEN=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('match_length',29))")
    MASH_KMER=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('mash_kmer',31))")
    PATH_JUMP=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('path_jump_max',0))")
    EDGE_JUMP=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('edge_jump_max',0))")
    PGGB_IMAGE=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('image','ghcr.io/pangenome/pggb:latest'))")
    CONFIG_SHA256=$(echo "$CFG_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('config_sha256',''))")

    echo "PGGB params (identical to chunk):"
    echo "  threads=$THREADS -p $MIN_ID -s $SEG_LEN -K $MASH_KMER -k $MATCH_LEN -j $PATH_JUMP -e $EDGE_JUMP"
    echo "  image=$PGGB_IMAGE"
    echo "  instance: ${PGGB_INSTANCE_TYPE:-mem3_ssd1_v2_x16}"

    docker pull "$PGGB_IMAGE" 2>&1 | tail -1

    INSTANCE_TYPE="${PGGB_INSTANCE_TYPE:-mem3_ssd1_v2_x16}"
    START_TS=$(date -u +"%Y-%m-%dT%H:%M:%S")
    START=$(date +%s)
    mkdir -p output

    docker run --rm \
        -v "$PWD/input.fa":/data/input.fa:ro \
        -v "$PWD/output":/data/output \
        "$PGGB_IMAGE" \
        pggb \
            -i /data/input.fa \
            -o /data/output \
            -t "$THREADS" \
            -n "$NUM_PATHS" \
            -p "$MIN_ID" \
            -s "$SEG_LEN" \
            -K "$MASH_KMER" \
            -k "$MATCH_LEN" \
            -j "$PATH_JUMP" \
            -e "$EDGE_JUMP" \
            2>&1 | tee pggb.log

    END=$(date +%s)
    END_TS=$(date -u +"%Y-%m-%dT%H:%M:%S")
    DURATION=$((END - START))
    echo "PGGB finished in ${DURATION}s"

    # Locate exactly one *final.gfa (FATAL if zero or >1)
    GFA_FILE=$(find output -name "*final.gfa" -type f | head -1)
    if [ -z "$GFA_FILE" ]; then
        GFA_FILE=$(find output -name "*final.gfa" -type f 2>/dev/null | head -1)
    fi
    if [ -z "$GFA_FILE" ]; then
        echo "FATAL: No *final.gfa produced"
        ls -la output/ 2>/dev/null || true
        exit 1
    fi
    GFA_COUNT=$(find output -name "*final.gfa" -type f 2>/dev/null | wc -l)
    if [ "$GFA_COUNT" -gt 1 ]; then
        echo "FATAL: Multiple *final.gfa found"
        exit 1
    fi

    cp "$GFA_FILE" baseline.gfa
    GFA_SIZE=$(wc -c < baseline.gfa)
    echo "  Final GFA: baseline.gfa ($GFA_SIZE bytes)"

    IMAGE_DIGEST=$(docker inspect "$PGGB_IMAGE" 2>/dev/null | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    digest = data[0]['RepoDigests'][0] if data[0].get('RepoDigests') else 'unknown'
    print(digest)
except: print('unknown')
" 2>/dev/null || echo "unknown")

    cat > metadata.json << JSONEOF
{
  "method": "pggb_baseline",
  "container": "$PGGB_IMAGE",
  "image_digest": "$IMAGE_DIGEST",
  "num_paths": $NUM_PATHS,
  "threads": $THREADS,
  "instance_type": "$INSTANCE_TYPE",
  "wall_seconds": $DURATION,
  "start_timestamp": "$START_TS",
  "stop_timestamp": "$END_TS",
  "final_gfa_size_bytes": $GFA_SIZE,
  "status": "completed",
  "config_sha256": "$CONFIG_SHA256",
  "pggb_params": {
    "minimum_identity": $MIN_ID,
    "segment_length": $SEG_LEN,
    "match_length": $MATCH_LEN,
    "mash_kmer": $MASH_KMER,
    "path_jump_max": $PATH_JUMP,
    "edge_jump_max": $EDGE_JUMP
  }
}
JSONEOF

    GFA_ID=$(dx upload baseline.gfa --brief)
    LOG_ID=$(dx upload pggb.log --brief)
    META_ID=$(dx upload metadata.json --brief)

    dx-jobutil-add-output gfa "$GFA_ID"
    dx-jobutil-add-output log "$LOG_ID"
    dx-jobutil-add-output metadata "$META_ID"

    echo "=== Done ==="
}