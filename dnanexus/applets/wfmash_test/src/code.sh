#!/usr/bin/env bash
set -euo pipefail
IMAGE="ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474"
docker pull "$IMAGE" 2>&1 | tail -1
dx download "$fasta" -o input.fa
mkdir -p out

echo "INPUT_SHA256=$(sha256sum input.fa | awk '{print $1}')"

echo "=== PASS1 ==="
docker run --rm -v "$PWD":/data "$IMAGE" bash -c '
wfmash -s 5000 -l 25000 -p 90 -n 1 -k 19 -H 0.001 -Y "#" -t 8 --tmp-base /data/out /data/input.fa --hg-filter-ani-diff 30 --approx-map > /data/out/mappings.paf 2> /data/out/p1.err
'
echo "PASS1=$(wc -l < out/mappings.paf) rec"

echo "=== MINIMAL ==="
docker run --rm -v "$PWD":/data "$IMAGE" bash -c '
wfmash -i /data/out/mappings.paf /data/input.fa > /data/out/minimal.paf 2> /data/out/min.err
'
echo "MINIMAL=$(wc -l < out/minimal.paf) rec, $(wc -c < out/minimal.paf) bytes"
echo "MINIMAL_CG=$(grep -c "cg:Z:" out/minimal.paf 2>/dev/null || echo 0)"
head -3 out/minimal.paf 2>/dev/null || echo "(empty)"
echo "=== MIN STDERR ==="
head -10 out/min.err 2>/dev/null || echo "(none)"

echo "=== CONTROL ==="
docker run --rm -v "$PWD":/data "$IMAGE" bash -c '
wfmash -s 5000 -l 25000 -p 90 -n 1 -k 19 -H 0.001 -Y "#" -t 8 --tmp-base /data/out /data/input.fa --lower-triangular --hg-filter-ani-diff 30 -i /data/out/mappings.paf --invert-filtering > /data/out/control.paf 2> /data/out/ctrl.err
'
echo "CONTROL=$(wc -l < out/control.paf) rec, $(wc -c < out/control.paf) bytes"
echo "CONTROL_CG=$(grep -c "cg:Z:" out/control.paf 2>/dev/null || echo 0)"

echo ""
echo "=== RESULTS ==="
echo "PASS1=$(wc -l < out/mappings.paf) rec"
echo "MINIMAL=$(wc -l < out/minimal.paf) rec, cg=$(grep -c 'cg:Z:' out/minimal.paf 2>/dev/null || echo 0)"
echo "CONTROL=$(wc -l < out/control.paf) rec, cg=$(grep -c 'cg:Z:' out/control.paf 2>/dev/null || echo 0)"
echo "DONE"

# Publish PAFs (optional — job succeeds even if empty)
for name in mappings minimal control; do
    local_file="out/${name}.paf"
    if [ -f "$local_file" ]; then
        fid=$(dx upload "$local_file" --brief 2>/dev/null || echo "")
        if [ -n "$fid" ]; then
            dx-jobutil-add-output "${name}_paf" "$fid" 2>/dev/null || true
            echo "PUBLISHED ${name}_paf"
        fi
    else
        echo "SKIP ${name}_paf (not found)"
    fi
done
