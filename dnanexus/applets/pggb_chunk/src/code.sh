#!/usr/bin/env bash
# pggb_chunk applet — runs PGGB via Docker inside a DNAnexus job.
# All PGGB parameters read from pggb_config_json (single source of truth).
# Image pinned to a specific digest — NEVER :latest.
set -e -o pipefail

main() {
    echo "=== PGGB Chunk Graph Builder ==="
    echo "Input FASTA: $fasta_name"
    echo ""

    # Download input file
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

    echo "PGGB params from pggb_config_json:"
    echo "  threads=$THREADS -p $MIN_ID -s $SEG_LEN -K $MASH_KMER -k $MATCH_LEN -j $PATH_JUMP -e $EDGE_JUMP"
    echo "  image=$PGGB_IMAGE"
    echo ""

    # Pull/ensure PGGB Docker image
    docker pull "$PGGB_IMAGE" 2>&1 | tail -1

    # Run PGGB
    INSTANCE_TYPE="${PGGB_INSTANCE_TYPE:-mem3_ssd1_v2_x16}"
    START_TS=$(date -u +"%Y-%m-%dT%H:%M:%S")
    START=$(date +%s)
    mkdir -p output

    # PGGB requires samtools faidx next to the FASTA. Mount the workdir so
    # the .fai can be written, then index + build in one container.
    docker run --rm \
        -v "$PWD":/data \
        "$PGGB_IMAGE" \
        bash -lc "samtools faidx /data/input.fa && pggb \
            -i /data/input.fa \
            -o /data/output \
            -t $THREADS \
            -n $NUM_PATHS \
            -p $MIN_ID \
            -s $SEG_LEN \
            -K $MASH_KMER \
            -k $MATCH_LEN \
            -j $PATH_JUMP \
            -e $EDGE_JUMP" \
            2>&1 | tee pggb.log

    END=$(date +%s)
    END_TS=$(date -u +"%Y-%m-%dT%H:%M:%S")
    DURATION=$((END - START))
    echo "PGGB finished in ${DURATION}s"

    # Locate exactly one *final.gfa (FATAL if zero or >1)
    GFA_FILE=$(find output -name "*final.gfa" -type f | head -1)
    if [ -z "$GFA_FILE" ]; then
        # PGGB may put it in a subdirectory
        GFA_FILE=$(find output -name "*final.gfa" -type f 2>/dev/null | head -1)
    fi
    if [ -z "$GFA_FILE" ]; then
        echo "FATAL: No *final.gfa produced by PGGB"
        ls -la output/ 2>/dev/null || true
        find output -name "*.gfa" -type f 2>/dev/null || true
        exit 1
    fi
    # Check for duplicates
    GFA_COUNT=$(find output -name "*final.gfa" -type f 2>/dev/null | wc -l)
    if [ "$GFA_COUNT" -gt 1 ]; then
        echo "FATAL: Multiple *final.gfa found ($GFA_COUNT)"
        find output -name "*final.gfa" -type f 2>/dev/null
        exit 1
    fi

    cp "$GFA_FILE" merged.gfa
    GFA_SIZE=$(wc -c < merged.gfa)
    echo "  Final GFA: merged.gfa ($GFA_SIZE bytes)"

    # Get image digest
    IMAGE_DIGEST=$(docker inspect "$PGGB_IMAGE" 2>/dev/null | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    digest = data[0]['RepoDigests'][0] if data[0].get('RepoDigests') else 'unknown'
    print(digest)
except: print('unknown')
" 2>/dev/null || echo "unknown")

    # Write proper metadata
    cat > metadata.json << JSONEOF
{
  "method": "pggb_chunk",
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
  },
  "full_command": "pggb -i /data/input.fa -o /data/output -t $THREADS -n $NUM_PATHS -p $MIN_ID -s $SEG_LEN -K $MASH_KMER -k $MATCH_LEN -j $PATH_JUMP -e $EDGE_JUMP"
}
JSONEOF

    # Upload results via job output links
    GFA_ID=$(dx upload merged.gfa --brief)
    LOG_ID=$(dx upload pggb.log --brief)
    META_ID=$(dx upload metadata.json --brief)

    dx-jobutil-add-output gfa "$GFA_ID"
    dx-jobutil-add-output log "$LOG_ID"
    dx-jobutil-add-output metadata "$META_ID"

    echo "=== Done ==="
}
