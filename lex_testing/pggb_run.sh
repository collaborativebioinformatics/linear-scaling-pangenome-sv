#!/usr/bin/env bash
# pggb_run.sh — Run PGGB on one FASTA and record what it cost.
#
# Measures wall-clock time and peak container memory, then reports the
# resulting graph's size. Appends one row to lex_testing/metrics/runs.tsv.
#
# Writes ONLY inside lex_testing/. Never touches DNAnexus.
#
# Usage:
#   bash lex_testing/pggb_run.sh <input.fa> <label> [threads]
#
# Example:
#   bash lex_testing/pggb_run.sh lex_testing/inputs/smoke_1mb.fa smoke_all5 16

set -euo pipefail

INPUT="${1:?Usage: pggb_run.sh <input.fa> <label> [threads]}"
LABEL="${2:?Need a label for this run}"
THREADS="${3:-16}"

PGGB_IMAGE="${PGGB_IMAGE:-ghcr.io/pangenome/pggb:latest}"
P_IDENT="${P_IDENT:-90}"     # -p  minimum alignment identity
S_SEG="${S_SEG:-5000}"       # -s  segment length
K_MIN="${K_MIN:-29}"         # -k  min match length

RUN_DIR="lex_testing/runs/${LABEL}"
METRICS_DIR="lex_testing/metrics"
mkdir -p "$RUN_DIR" "$METRICS_DIR"
TSV="$METRICS_DIR/runs.tsv"

[ -f "$INPUT" ] || { echo "ERROR: no such input: $INPUT"; exit 1; }

# ── Prepare a bgzipped, indexed copy (in our own inputs dir, never DNAnexus) ──
if [[ "$INPUT" == *.gz ]]; then
    GZ="$INPUT"
else
    GZ="${INPUT}.gz"
    [ -f "$GZ" ] || { echo "bgzip -> $GZ"; bgzip -@ "$THREADS" -k "$INPUT"; }
fi
[ -f "${GZ}.fai" ] && [ -f "${GZ}.gzi" ] || samtools faidx "$GZ"

NSEQ=$(wc -l < "${GZ}.fai" | tr -d ' ')
INBP=$(awk '{s+=$2} END{print s+0}' "${GZ}.fai")

echo "=============================================="
echo " PGGB run: $LABEL"
echo " input     : $GZ"
echo " sequences : $NSEQ   (${INBP} bp total)"
echo " threads   : $THREADS"
echo " params    : -p $P_IDENT -s $S_SEG -k $K_MIN"
echo "=============================================="

docker image inspect "$PGGB_IMAGE" &>/dev/null || docker pull "$PGGB_IMAGE"

IN_DIR="$(cd "$(dirname "$GZ")" && pwd)"
OUT_DIR="$(cd "$RUN_DIR" && pwd)"
LOG="$RUN_DIR/pggb.log"
CIDFILE="$RUN_DIR/.cid"; rm -f "$CIDFILE"

# ── Launch, then poll memory while it runs ──────────────────────────────────
START=$(date +%s)
docker run --rm --cidfile "$CIDFILE" \
    -v "$IN_DIR":/data/input:ro \
    -v "$OUT_DIR":/data/output \
    "$PGGB_IMAGE" \
    pggb -i "/data/input/$(basename "$GZ")" \
         -o /data/output \
         -t "$THREADS" -n "$NSEQ" \
         -p "$P_IDENT" -s "$S_SEG" -k "$K_MIN" \
         -j 0 -e 0 > "$LOG" 2>&1 &
DPID=$!

to_mib() {
    awk '{v=$1
          if (v ~ /GiB$/)      {sub(/GiB$/,"",v); v*=1024}
          else if (v ~ /MiB$/) {sub(/MiB$/,"",v)}
          else if (v ~ /KiB$/) {sub(/KiB$/,"",v); v/=1024}
          else if (v ~ /B$/)   {sub(/B$/,"",v);   v/=1048576}
          printf "%.1f", v}'
}

PEAK=0
while kill -0 "$DPID" 2>/dev/null; do
    if [ -s "$CIDFILE" ]; then
        CUR=$(docker stats --no-stream --format '{{.MemUsage}}' \
                "$(cat "$CIDFILE")" 2>/dev/null | to_mib || true)
        if [ -n "${CUR:-}" ]; then
            PEAK=$(awk -v a="$PEAK" -v b="$CUR" 'BEGIN{print (b>a)?b:a}')
        fi
    fi
    sleep 5
done
RC=0; wait "$DPID" || RC=$?
END=$(date +%s)
WALL=$((END - START))
rm -f "$CIDFILE"

# ── Locate the final graph (prefer the smoothed final GFA) ──────────────────
GFA=""
for PAT in "*.smooth.final.gfa" "*final*.gfa" "*.gfa"; do
    GFA=$(find "$RUN_DIR" -name "$PAT" -type f 2>/dev/null | head -1)
    [ -n "$GFA" ] && break
done

if [ -z "$GFA" ]; then
    echo "FAILED after ${WALL}s (exit $RC) — no GFA produced. See $LOG"
    STATUS="FAILED"; NODES=0; EDGES=0; GBP=0; GFA="-"
else
    cp "$GFA" "$RUN_DIR/${LABEL}.gfa"
    GFA="$RUN_DIR/${LABEL}.gfa"
    STATUS="OK"
    NODES=$(grep -c '^S' "$GFA" || true)
    EDGES=$(grep -c '^L' "$GFA" || true)
    GBP=$(awk -F'\t' '$1=="S"{s+=length($3)} END{print s+0}' "$GFA")
fi

[ -f "$TSV" ] || printf 'label\tsequences\tinput_bp\tthreads\twall_s\tpeak_mib\tnodes\tedges\tgraph_bp\tstatus\ttimestamp\n' > "$TSV"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$LABEL" "$NSEQ" "$INBP" "$THREADS" "$WALL" "$PEAK" \
    "$NODES" "$EDGES" "$GBP" "$STATUS" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$TSV"

echo "----------------------------------------------"
printf ' %-12s %s\n' "status"  "$STATUS"
printf ' %-12s %ss (%s)\n' "wall time" "$WALL" "$(printf '%dh%02dm%02ds' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60)))"
printf ' %-12s %s MiB\n' "peak mem" "$PEAK"
printf ' %-12s %s nodes / %s edges / %s bp\n' "graph" "$NODES" "$EDGES" "$GBP"
printf ' %-12s %s\n' "gfa" "$GFA"
printf ' %-12s %s\n' "metrics" "$TSV"
echo "----------------------------------------------"
exit "$RC"
