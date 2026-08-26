#!/usr/bin/env bash
# benchmark.sh — Measure how PGGB's cost scales with the number of haplotypes.
#
# Three modes:
#   pairwise    reference + ONE haplotype, once per haplotype.
#               N separate small graphs. Cost grows LINEARLY with samples.
#   cumulative  reference + 1, +2, +3, ... haplotypes.
#               Traces the growth curve as the input widens.
#   all         every sequence in one graph (the monolithic baseline).
#               All-vs-all alignment: cost grows QUADRATICALLY.
#
# Writes ONLY inside lex_testing/. Never touches DNAnexus.
#
# Usage:
#   bash lex_testing/benchmark.sh <input.fa> <mode> [threads] [tag]
#   bash lex_testing/benchmark.sh lex_testing/inputs/smoke_1mb.fa cumulative 16 smoke

set -euo pipefail

INPUT="${1:?Usage: benchmark.sh <input.fa> <pairwise|cumulative|all> [threads] [tag]}"
MODE="${2:?Need a mode: pairwise | cumulative | all}"
THREADS="${3:-16}"
TAG="${4:-$(basename "${INPUT%.fa}")}"

SUB_DIR="lex_testing/inputs/subsets"
mkdir -p "$SUB_DIR"
[ -f "${INPUT}.fai" ] || samtools faidx "$INPUT"

mapfile -t SEQS < <(cut -f1 "${INPUT}.fai")
[ "${#SEQS[@]}" -ge 2 ] || { echo "ERROR: need >=2 sequences"; exit 1; }

REF="${SEQS[0]}"
HAPS=("${SEQS[@]:1}")
echo "Reference   : $REF"
echo "Haplotypes  : ${#HAPS[@]}"
printf '  %s\n' "${HAPS[@]}"
echo

# short, filesystem-safe name from a PanSN header (HG00673#2#JAH... -> HG00673_2)
short() { echo "$1" | awk -F'#' '{print (NF>=2)? $1"_"$2 : $1}' | tr -c 'A-Za-z0-9_.-' '_'; }

make_subset() {   # $1 = output fasta, rest = sequence names
    local out="$1"; shift
    if [ -f "$out" ]; then echo "  reuse $(basename "$out")"; return; fi
    samtools faidx "$INPUT" "$@" > "$out"
    samtools faidx "$out"
}

case "$MODE" in
  pairwise)
    for h in "${HAPS[@]}"; do
        lbl="${TAG}_pair_$(short "$h")"
        f="$SUB_DIR/${lbl}.fa"
        echo "--- $lbl ---"
        make_subset "$f" "$REF" "$h"
        bash lex_testing/pggb_run.sh "$f" "$lbl" "$THREADS" || echo "  (continuing)"
    done
    ;;
  cumulative)
    for ((k=1; k<=${#HAPS[@]}; k++)); do
        lbl="${TAG}_cum_${k}hap"
        f="$SUB_DIR/${lbl}.fa"
        echo "--- $lbl (reference + $k haplotype(s)) ---"
        make_subset "$f" "$REF" "${HAPS[@]:0:$k}"
        bash lex_testing/pggb_run.sh "$f" "$lbl" "$THREADS" || echo "  (continuing)"
    done
    ;;
  all)
    echo "--- ${TAG}_all (${#SEQS[@]} sequences together) ---"
    bash lex_testing/pggb_run.sh "$INPUT" "${TAG}_all" "$THREADS"
    ;;
  *) echo "Unknown mode: $MODE"; exit 1 ;;
esac

echo
bash lex_testing/summary.sh
