#!/usr/bin/env bash
# fetch_inputs.sh — Pull benchmark inputs from DNAnexus.
#
# READ-ONLY. This script only ever calls `dx download`. It never uploads,
# creates folders, or deletes anything in the DNAnexus project.
#
# Usage: bash lex_testing/fetch_inputs.sh [smoke|full|both]

set -euo pipefail

WHICH="${1:-smoke}"
DEST="lex_testing/inputs"
mkdir -p "$DEST"

: "${DXPROJ:?Set DXPROJ first, e.g. export DXPROJ=project-JB6zQBj0ZQv2Bk79ggBBv76Z}"

SMOKE_ID="file-JB7Zkz00ZQvJKQ7j955zxJQJ"   # chr21_20000000_21000000_multi.fa
FULL_ID="file-JB7YBYj0ZQv0yQZGqZkPx45G"    # chr21_multi.fa

get() {
    local fid="$1" out="$2"
    if [ -f "$out" ]; then
        echo "  EXISTS $out ($(du -h "$out" | cut -f1))"
        return
    fi
    echo "  Downloading $fid -> $out"
    dx download "${DXPROJ}:${fid}" -o "$out"
}

echo "=== Fetching benchmark inputs (read-only) ==="
case "$WHICH" in
    smoke) get "$SMOKE_ID" "$DEST/smoke_1mb.fa" ;;
    full)  get "$FULL_ID"  "$DEST/chr21_full.fa" ;;
    both)  get "$SMOKE_ID" "$DEST/smoke_1mb.fa"
           get "$FULL_ID"  "$DEST/chr21_full.fa" ;;
    *) echo "Usage: $0 [smoke|full|both]"; exit 1 ;;
esac

for f in "$DEST"/*.fa; do
    [ -f "$f" ] || continue
    [ -f "${f}.fai" ] || samtools faidx "$f"
    echo "--- $(basename "$f") ---"
    awk '{printf "    %-32s %12d bp\n", $1, $2}' "${f}.fai"
done
echo "=== Done. Nothing on DNAnexus was modified. ==="
