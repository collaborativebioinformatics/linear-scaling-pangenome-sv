"""Write a diagnostic applet that runs only wfmash test matrix, uploads all PAFs."""
import os

script = """#!/usr/bin/env bash
set -euo pipefail

IMAGE="ghcr.io/pangenome/pggb@sha256:44f3563c5c3e2032c293fe309eed178b754021b07963ca4857e4ae5efda86474"
docker pull "$IMAGE" 2>&1 | tail -1

dx download "$fasta" -o input.fa
samtools faidx input.fa

mkdir -p /tmp/out
for mode in A B C D E1 E2; do mkdir -p "/tmp/out/$mode"; done

# A: simple all-vs-all, no lower-triangular
echo "=== MODE A ==="
docker run --rm -v /tmp/out/A:/out -v "$PWD/input.fa":/in.fa:ro "$IMAGE" bash -c '
wfmash -s 5000 -l 25000 -p 90 -n 1 -k 19 -H 0.001 -Y "#" -t 8 --tmp-base /tmp /in.fa --hg-filter-ani-diff 30 > /out/A.paf 2> /out/A.stderr
'
wc -c /tmp/out/A/A.paf; wc -l /tmp/out/A/A.paf

# B: with --lower-triangular (V5 mode)
echo "=== MODE B ==="
docker run --rm -v /tmp/out/B:/out -v "$PWD/input.fa":/in.fa:ro "$IMAGE" bash -c '
wfmash -s 5000 -l 25000 -p 90 -n 1 -k 19 -H 0.001 -Y "#" -t 8 --tmp-base /tmp /in.fa --lower-triangular --hg-filter-ani-diff 30 > /out/B.paf 2> /out/B.stderr
'
wc -c /tmp/out/B/B.paf; wc -l /tmp/out/B/B.paf

# C: no group filter
echo "=== MODE C ==="
docker run --rm -v /tmp/out/C:/out -v "$PWD/input.fa":/in.fa:ro "$IMAGE" bash -c '
wfmash -s 5000 -l 25000 -p 90 -n 1 -k 19 -H 0.001 -t 8 --tmp-base /tmp /in.fa --hg-filter-ani-diff 30 > /out/C.paf 2> /out/C.stderr
'
wc -c /tmp/out/C/C.paf; wc -l /tmp/out/C/C.paf

# D: no group filter + lower-triangular
echo "=== MODE D ==="
docker run --rm -v /tmp/out/D:/out -v "$PWD/input.fa":/in.fa:ro "$IMAGE" bash -c '
wfmash -s 5000 -l 25000 -p 90 -n 1 -k 19 -H 0.001 -t 8 --tmp-base /tmp /in.fa --lower-triangular --hg-filter-ani-diff 30 > /out/D.paf 2> /out/D.stderr
'
wc -c /tmp/out/D/D.paf; wc -l /tmp/out/D/D.paf

# E1: PGGB pass1 --approx-map
echo "=== MODE E1 ==="
docker run --rm -v /tmp/out/E1:/out -v "$PWD/input.fa":/in.fa:ro "$IMAGE" bash -c '
wfmash -s 5000 -l 25000 -p 90 -n 1 -k 19 -H 0.001 -Y "#" -t 8 --tmp-base /tmp /in.fa --hg-filter-ani-diff 30 --approx-map > /out/mappings.paf 2> /out/E1.stderr
'
wc -c /tmp/out/E1/mappings.paf; wc -l /tmp/out/E1/mappings.paf

# E2: PGGB pass2 --lower-triangular --invert-filtering -i mappings.paf
echo "=== MODE E2 ==="
docker run --rm -v /tmp/out/E2:/out -v /tmp/out/E1/mappings.paf:/mappings.paf:ro -v "$PWD/input.fa":/in.fa:ro "$IMAGE" bash -c '
wfmash -s 5000 -l 25000 -p 90 -n 1 -k 19 -H 0.001 -Y "#" -t 8 --tmp-base /tmp /in.fa --lower-triangular --hg-filter-ani-diff 30 -i /mappings.paf --invert-filtering > /out/alignments.paf 2> /out/E2.stderr
'
wc -c /tmp/out/E2/alignments.paf; wc -l /tmp/out/E2/alignments.paf

# Summary
echo ""
echo "============================================"
echo "RESULTS"
echo "============================================"
for mode in A B C D; do
    f="/tmp/out/${mode}/${mode}.paf"
    echo "$mode: $(wc -l < $f) rec  $(wc -c < $f) bytes"
done
echo "E1 (mapping): $(wc -l < /tmp/out/E1/mappings.paf) rec"
echo "E2 (detailed): $(wc -l < /tmp/out/E2/alignments.paf) rec"

# Upload all as outputs
A_ID=$(dx upload /tmp/out/A/A.paf --brief)
B_ID=$(dx upload /tmp/out/B/B.paf --brief)
C_ID=$(dx upload /tmp/out/C/C.paf --brief)
D_ID=$(dx upload /tmp/out/D/D.paf --brief)
E1_ID=$(dx upload /tmp/out/E1/mappings.paf --brief)
E2_ID=$(dx upload /tmp/out/E2/alignments.paf --brief)

dx-jobutil-add-output result "$A_ID"
"""

os.makedirs("dnanexus/applets/wfmash_diag/src", exist_ok=True)
with open("dnanexus/applets/wfmash_diag/src/code.sh", "w") as f:
    f.write(script)

dxapp = {
    "name": "wfmash_diag", "title": "WFMASH Diagnostic",
    "version": "0.1.0",
    "inputSpec": [{"name": "fasta", "class": "file", "help": "Chunk FASTA"}],
    "outputSpec": [{"name": "result", "class": "file", "help": "Results"}],
    "runSpec": {"interpreter": "bash", "file": "src/code.sh", "distribution": "Ubuntu", "release": "24.04",
                "systemRequirements": {"*": {"instanceType": "mem3_ssd1_v2_x16"}}},
    "access": {"network": ["*"]}
}
with open("dnanexus/applets/wfmash_diag/dxapp.json", "w") as f:
    import json; json.dump(dxapp, f, indent=2)
print("Applet dir written")