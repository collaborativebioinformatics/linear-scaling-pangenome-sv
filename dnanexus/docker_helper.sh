#!/usr/bin/env bash
# docker_helper.sh — Run bioinformatics tools via Docker containers.
# Source this file from other scripts: source dnanexus/docker_helper.sh
#
# Usage:
#   source dnanexus/docker_helper.sh
#   run_pggb input.fa output_dir threads num_paths
#   run_vg <vg args>

set -uo pipefail

PGGB_IMAGE="${PGGB_IMAGE:-ghcr.io/pangenome/pggb:latest}"
VG_IMAGE="${VG_IMAGE:-quay.io/vgteam/vg:v1.74.1}"

ensure_image() {
    local image="$1"
    if ! docker image inspect "$image" &>/dev/null; then
        echo "Pulling Docker image: $image"
        docker pull "$image"
    fi
}

run_pggb() {
    local input="$1"
    local outdir="$2"
    local threads="${3:-16}"
    local num_paths="${4:-5}"

    ensure_image "$PGGB_IMAGE"

    local input_dir; input_dir="$(cd "$(dirname "$input")" && pwd)"
    local input_file; input_file="$(basename "$input")"
    local outdir_abs; outdir_abs="$(cd "$(dirname "$outdir")" && pwd)/$(basename "$outdir")"
    mkdir -p "$outdir"

    echo "Running PGGB via Docker..."
    echo "  Image: $PGGB_IMAGE"
    echo "  Paths: $num_paths"

    docker run --rm \
        -v "$input_dir":/data/input:ro \
        -v "$outdir_abs":/data/output \
        "$PGGB_IMAGE" \
        pggb \
            -i "/data/input/$input_file" \
            -o "/data/output" \
            -t "$threads" \
            -n "$num_paths" \
            -p 90 \
            -s 5000 \
            -k 29 \
            -w 50000 \
            -j 0 \
            -e 0

    # Find the output GFA
    local gfa; gfa=$(find "$outdir" -name "*.gfa" -type f | head -1)
    if [ -z "$gfa" ]; then
        echo "WARNING: PGGB produced no GFA in $outdir"
        return 1
    fi
    echo "$gfa"
}

run_vg() {
    ensure_image "$VG_IMAGE"
    docker run --rm -v "$PWD":/data "$VG_IMAGE" vg "$@"
}

run_odgi() {
    ensure_image "$PGGB_IMAGE"  # odgi is bundled inside pggb container
    docker run --rm -v "$PWD":/data "$PGGB_IMAGE" odgi "$@"
}

echo "Docker helper loaded."
echo "  PGGB: $PGGB_IMAGE"
echo "  vg:   $VG_IMAGE"