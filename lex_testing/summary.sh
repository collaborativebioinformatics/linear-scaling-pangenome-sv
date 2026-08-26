#!/usr/bin/env bash
# summary.sh — Print the cost table collected so far.
set -euo pipefail
TSV="lex_testing/metrics/runs.tsv"
[ -f "$TSV" ] || { echo "No runs yet ($TSV missing)."; exit 0; }

echo "=== PGGB cost table ==="
column -t -s "$(printf '\t')" "$TSV"
echo
echo "=== Cost per haplotype (wall seconds / sequences) ==="
awk -F'\t' 'NR==1{next} $10=="OK" && $2>0 {
    printf "  %-28s %2d seq  %7ds  %9.1f s/seq  %8s MiB\n", $1, $2, $5, $5/$2, $6
}' "$TSV"
