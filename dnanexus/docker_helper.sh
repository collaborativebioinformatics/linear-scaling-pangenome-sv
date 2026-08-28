#!/usr/bin/env bash
# docker_helper.sh — Run bioinformatics tools via Docker containers.
# Source this file from other scripts: source dnanexus/docker_helper.sh
#
# LEGACY/DEVELOPMENT_ONLY helper for LOCAL execution.
# The canonical PGGB execution path is the DNAnexus applets
# (dnanexus/applets/pggb_chunk + pggb_baseline), which read pggb_config_json.
# This helper mirrors the SAME canonical parameters from gen_pggb_config.py.

set -uo pipefail

read_pggb_config() {
    # Read the canonical PGGB image from the single source of truth.
    python3 scripts/gen_pggb_config.py 2>/dev/null | \
        python3 -c "import sys,json; print(json.load(sys.stdin)['image'])" 2>/dev/null \
        || python3 -c "import yaml; print(yaml.safe_load(open('config/pipeline.yaml'))['pggb']['image'])" 2>/dev/null \
        || echo "ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474"
}
PGGB_IMAGE="${PGGB_IMAGE:-$(read_pggb_config)}"
VG_IMAGE="${VG_IMAGE:-quay.io/vgteam/vg:v1.74.1}"

ensure_image() {
    local image="$1"
    if ! docker image inspect "$image" &>/dev/null; then
        echo "Pulling Docker image: $image"
        docker pull "$image"
    fi
}

# run_pggb: local PGGB run. Parameters come from the canonical config JSON,
# never from hard-coded flags. This mirrors the applet code.sh behavior.
run_pggb() {
    local input="$1"
    local outdir="$2"
    local threads="${3:-}"
    local num_paths="${4:-}"

    ensure_image "$PGGB_IMAGE"

    local input_dir; input_dir="$(cd "$(dirname "$input")" && pwd)"
    local input_file; input_file="$(basename "$input")"
    local outdir_abs; outdir_abs="$(cd "$(dirname "$outdir")" && pwd)/$(basename "$outdir")"
    mkdir -p "$outdir"

    # Load canonical params (defaults match config/pipeline.yaml)
    local cfg
    cfg="$(python3 scripts/gen_pggb_config.py 2>/dev/null || echo '{}')"
    local MIN_ID; MIN_ID="$(echo "$cfg" | python3 -c "import sys,json;print(json.load(sys.stdin).get('minimum_identity',90))" 2>/dev/null || echo 90)"
    local SEG_LEN; SEG_LEN="$(echo "$cfg" | python3 -c "import sys,json;print(json.load(sys.stdin).get('segment_length',5000))" 2>/dev/null || echo 5000)"
    local MASH_KMER; MASH_KMER="$(echo "$cfg" | python3 -c "import sys,json;print(json.load(sys.stdin).get('mash_kmer',31))" 2>/dev/null || echo 31)"
    local MATCH_LEN; MATCH_LEN="$(echo "$cfg" | python3 -c "import sys,json;print(json.load(sys.stdin).get('match_length',29))" 2>/dev/null || echo 29)"
    local PATH_JUMP; PATH_JUMP="$(echo "$cfg" | python3 -c "import sys,json;print(json.load(sys.stdin).get('path_jump_max',0))" 2>/dev/null || echo 0)"
    local EDGE_JUMP; EDGE_JUMP="$(echo "$cfg" | python3 -c "import sys,json;print(json.load(sys.stdin).get('edge_jump_max',0))" 2>/dev/null || echo 0)"
    local THREADS; THREADS="${threads:-$(echo "$cfg" | python3 -c "import sys,json;print(json.load(sys.stdin).get('threads',8))" 2>/dev/null || echo 8)}"
    local NP; NP="${num_paths:-$(grep -c '^>' "$input")}"

    echo "Running PGGB via Docker..."
    echo "  Image: $PGGB_IMAGE"
    echo "  Paths: $NP  Threads: $THREADS"

    docker run --rm \
        -v "$input_dir":/data/input:ro \
        -v "$outdir_abs":/data/output \
        "$PGGB_IMAGE" \
        bash -lc "samtools faidx /data/input/$input_file && pggb \
            -i /data/input/$input_file \
            -o /data/output \
            -t $THREADS \
            -n $NP \
            -p $MIN_ID \
            -s $SEG_LEN \
            -K $MASH_KMER \
            -k $MATCH_LEN \
            -j $PATH_JUMP \
            -e $EDGE_JUMP"

    # Find the output GFA
    local gfa; gfa=$(find "$outdir" -name "*final.gfa" -type f 2>/dev/null | head -1)
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