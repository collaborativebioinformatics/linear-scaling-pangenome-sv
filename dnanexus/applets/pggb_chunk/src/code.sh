#!/usr/bin/env bash
# pggb_chunk applet — runs PGGB via Docker inside a DNAnexus job.
set -e -o pipefail

main() {
    echo "=== PGGB Chunk Graph Builder ==="
    echo "Input FASTA: $fasta_name"
    echo ""

    # Download input file
    dx download "$fasta" -o input.fa
    NUM_PATHS=$(grep -c '^>' input.fa)
    echo "Paths in input: $NUM_PATHS"

    # Pull/ensure PGGB Docker image
    PGGB_IMAGE="ghcr.io/pangenome/pggb:latest"
    docker pull "$PGGB_IMAGE" 2>&1 | tail -1 || true

    # Run PGGB
    START=$(date +%s)
    mkdir -p output

    docker run --rm \
        -v "$PWD/input.fa":/data/input.fa:ro \
        -v "$PWD/output":/data/output \
        "$PGGB_IMAGE" \
        pggb \
            -i /data/input.fa \
            -o /data/output \
            -t "$(nproc)" \
            -n "$NUM_PATHS" \
            -p 90 \
            -s 5000 \
            -k 29 \
            -j 0 \
            -e 0 \
            2>&1 | tee pggb.log

    END=$(date +%s)
    DURATION=$((END - START))
    echo "PGGB finished in ${DURATION}s"

    # Locate the GFA output
    GFA_FILE=$(find output -name "*.gfa" -type f | head -1)
    if [ -z "$GFA_FILE" ]; then
        echo "ERROR: No GFA produced by PGGB"
        ls -la output/ 2>/dev/null || true
        exit 1
    fi

    cp "$GFA_FILE" merged.gfa

    # Write metadata
    cat > metadata.json << JSONEOF
{
  "method": "pggb_chunk",
  "container": "$PGGB_IMAGE",
  "num_paths": $NUM_PATHS,
  "wall_seconds": $DURATION,
  "status": "completed"
}
JSONEOF

    # Upload results
    GFA_ID=$(dx upload merged.gfa --brief)
    LOG_ID=$(dx upload pggb.log --brief)
    META_ID=$(dx upload metadata.json --brief)

    dx-jobutil-add-output gfa "$GFA_ID"
    dx-jobutil-add-output log "$LOG_ID"
    dx-jobutil-add-output metadata "$META_ID"

    echo "=== Done ==="
}