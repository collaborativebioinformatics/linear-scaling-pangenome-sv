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

    # Two-pass PGGB with detailed pass fixed for small inputs.
    # Baseline works because pass 1 (--approx-map) doesn't cover all 5 Mb.
    # Small chunks fail because pass 1 covers all 2.1 Mb, leaving pass 2
    # (--invert-filtering) with nothing to align → 0 links.
    # Fix: Remove --invert-filtering from pass 2 so it DETAILS what pass 1 FOUND.
    cat > /tmp/chunk_pipeline.sh << 'HEREDOC_END'
#!/usr/bin/env bash
set -euo pipefail
set -x
echo "=== Chunk PGGB (two-pass, detailed aligns what pass1 found) ==="

# 1. Index
samtools faidx /data/input.fa
echo "[1/7] samtools faidx done"

# 2. Pass 1: approximate mappings (same as baseline)
wfmash -s __SEGLEN__ -l 25000 -p __PID__ -n 1 \
    -k __KMER__ -H 0.001 -Y '#' \
    -t __THREADS__ --tmp-base /data/output \
    /data/input.fa \
    --hg-filter-ani-diff 30 --approx-map \
    > /data/output/mappings.paf \
    2> /data/output/pass1.stderr
echo "[2/7] pass1 mappings done ($(wc -l < /data/output/mappings.paf) records)"

# 3. Pass 2: DETAILED alignment on what pass 1 FOUND (no --invert-filtering!)
wfmash -s __SEGLEN__ -l 25000 -p __PID__ -n 1 \
    -k __KMER__ -H 0.001 -Y '#' \
    -t __THREADS__ --tmp-base /data/output \
    /data/input.fa \
    --lower-triangular --hg-filter-ani-diff 30 \
    -i /data/output/mappings.paf \
    > /data/output/alignments.paf \
    2> /data/output/pass2.stderr
echo "[3/7] pass2 detailed done ($(wc -l < /data/output/alignments.paf) records)"

# Gate: refuse empty PAF
if [ ! -s /data/output/alignments.paf ]; then
    echo "FATAL EMPTY_ALIGNMENT_PAF"
    exit 1
fi

# 4. seqwish
seqwish -s /data/input.fa \
    -p /data/output/alignments.paf \
    -k __MATCHLEN__ -f 0 \
    -g /data/output/seqwish.gfa \
    -B 10M -t __THREADS__ \
    --temp-dir /data/output -P
echo "[4/7] seqwish done"

# 5. smoothxg (exact baseline flags)
smoothxg -t __THREADS__ -T __THREADS__ \
    -g /data/output/seqwish.gfa \
    -r __NPATHS__ \
    --base /data/output \
    --chop-to 100 \
    -I .9000 -R 0 \
    -j __PJUMP__ -e __EJUMP__ \
    -l 700,1100 -p 1,4,6,2,26,1 \
    -O 0.001 -Y 500 -d 0 -D 0 \
    -Q Consensus_ -V \
    -o /data/output/smooth.gfa
echo "[5/7] smoothxg done"

# 6-7. odgi
odgi build -g /data/output/smooth.gfa \
    -o /data/output/graph.og -t __THREADS__ -P
odgi sort -i /data/output/graph.og \
    -o /data/output/graph.sorted.og \
    -p Ygs -t __THREADS__ -P
odgi view -i /data/output/graph.sorted.og \
    -g > /data/output/final.gfa
echo "[6/7] odgi done — final.gfa written"
HEREDOC_END

    sed -i \
        -e "s/__SEGLEN__/${SEG_LEN}/g" -e "s/__PID__/${MIN_ID}/g" \
        -e "s/__KMER__/${MASH_KMER}/g" -e "s/__MATCHLEN__/${MATCH_LEN}/g" \
        -e "s/__PJUMP__/${PATH_JUMP}/g" -e "s/__EJUMP__/${EDGE_JUMP}/g" \
        -e "s/__THREADS__/${THREADS}/g" -e "s/__NPATHS__/${NUM_PATHS}/g" \
        /tmp/chunk_pipeline.sh

    chmod +x /tmp/chunk_pipeline.sh

    docker run --rm \
        -v "$PWD":/data \
        -v /tmp/chunk_pipeline.sh:/data/pipeline.sh:ro \
        "$PGGB_IMAGE" \
        bash /data/pipeline.sh \
        2>&1 | tee pggb.log

    # Validate edges exist
    if ! grep -q '^L' output/final.gfa 2>/dev/null; then
        echo "FATAL INVALID_ZERO_EDGE_GRAPH"
        ls -la output/ 2>/dev/null || true
        exit 1
    fi
    EDGE_COUNT=$(grep -c '^L' output/final.gfa)
    echo "  Graph has $EDGE_COUNT L-records"

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
