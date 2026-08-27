#!/usr/bin/env bash
# run_pipeline.sh — Run the full chr21 smoke-test pipeline.
# Auto-detects $HOME/.pipeline-env/bin venv if present.
# Flow: environment -> HPRC manifest -> stage inputs -> interval map ->
#       baseline PGGB (DNAnexus applet) -> chunk FASTAs ->
#       DNAnexus parallel chunk PGGB -> merge (diagnostic) ->
#       benchmark -> web JSON
set -euo pipefail

# Auto-detect pipeline virtualenv
if [ -d "$HOME/.pipeline-env/bin" ]; then
    export PATH="$HOME/.pipeline-env/bin:$PATH"
    echo "  Using $HOME/.pipeline-env/bin"
fi
if [ -d "$HOME/.local/bin" ]; then
    export PATH="$HOME/.local/bin:$PATH"
fi

UPLOAD="${1:-}"
echo "Pipeline starting (target: chr21:20000000-21000000, 0-based half-open)"
cd "$(cd "$(dirname "$0")/.." && pwd)"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="results/logs"; mkdir -p "$LOG_DIR"
exec > >(tee "$LOG_DIR/pipeline_${TIMESTAMP}.log") 2>&1

echo "[1/9] Environment check"
bash scripts/check_environment.sh || { echo "Failed"; exit 1; }

echo "[2/9] Fetching HPRC index"
python3 scripts/fetch_hprc_index.py

echo "[3/9] Staging inputs from DNAnexus"
bash dnanexus/stage_inputs.sh

echo "[3b/9] Missing assemblies -> HPRC S3 fallback"
MISSING=0
for F in HG00673_mat_hprc_r2_v1.0.1.fa.gz HG00673_pat_hprc_r2_v1.0.1.fa.gz \
         HG00733_mat_hprc_r2_v1.0.1.fa.gz HG00733_pat_hprc_r2_v1.0.1.fa.gz; do
    [ ! -f "work/downloads/$F" ] && MISSING=$((MISSING + 1))
done
if [ "$MISSING" -gt 0 ]; then
    python3 scripts/download_hprc.py --execute || true
    bash dnanexus/upload_inputs.sh
else
    echo "  All 4 present. No download needed."
fi

echo "[3c/9] Preparing reference"
bash scripts/prepare_reference.sh

echo "[4/9] Mapping and preparing chr21 sequences"
python3 pipeline/prepare/map_chromosome.py

echo "  Validating mapping report..."
MAPPING="results/preparation/sequence_mapping.tsv"
[ ! -f "$MAPPING" ] && { echo "FATAL: no mapping report"; exit 1; }
python3 -c "
import csv
rows = list(csv.DictReader(open('$MAPPING'), delimiter='\\t'))
mapped = [r for r in rows if r.get('status') == 'mapped' and r.get('source_start') and r.get('source_end')]
if len(mapped) != 4:
    print(f'FATAL: Expected 4 mapped with source coords, got {len(mapped)}')
    exit(1)
print('  All 4 haplotypes have valid source coordinates. Proceeding.')
"

python3 pipeline/prepare/prepare_sequences.py

INPUT_FASTA="results/preparation/chr21_multi.fa"
[ ! -f "$INPUT_FASTA" ] && { echo "FATAL: no chr21 multi-FASTA"; exit 1; }
NP=$(grep -c '^>' "$INPUT_FASTA")
[ "$NP" -ne 5 ] && { echo "FATAL: expected 5 paths, got $NP"; exit 1; }
echo "  $NP paths verified"

REF_LEN=$(python3 -c "
s = []
with open('$INPUT_FASTA') as f:
    for line in f:
        if line.startswith('>GRCh38'): continue
        if line.startswith('>'): break
        s.append(line.strip())
print(len(''.join(s)))
")
echo "  Reference sequence length: $REF_LEN bp"
if [ "$REF_LEN" -gt 2000000 ]; then
    echo "FATAL: Reference is $REF_LEN bp (expected ~1,000,000). Full chr21 not sliced."
    exit 1
fi

echo "[5/9] Baseline PGGB graph (DNAnexus applet)"
BASELINE_GFA="results/baseline/baseline.gfa"
mkdir -p results/baseline
if [ -n "${DX_PROJECT_ID:-${DX_PROJECT_CONTEXT_ID:-}}" ]; then
    echo "  Building pggb_baseline applet..."
    BASELINE_APPLET=$(cd dnanexus/applets/pggb_baseline && dx build --destination /applets/pggb_baseline/ --brief 2>/dev/null || echo "")
    if [ -n "$BASELINE_APPLET" ]; then
        echo "  Uploading input FASTA..."
        FASTA_ID=$(dx upload "$INPUT_FASTA" --destination /data/prepared/chr21_multi.fa --brief 2>/dev/null || echo "")
        if [ -n "$FASTA_ID" ]; then
            echo "  Launching baseline PGGB job..."
            BASELINE_JOB=$(dx run "$BASELINE_APPLET" \
                -i fasta="$FASTA_ID" \
                --instance-type "${DX_INSTANCE_TYPE:-mem3_ssd1_v2_x16}" \
                --name "PGGB-Baseline" \
                --destination /graphs/baseline/ \
                --brief 2>/dev/null || echo "")
            if [ -n "$BASELINE_JOB" ]; then
                echo "  Baseline job: $BASELINE_JOB"
                dx wait "$BASELINE_JOB" 2>/dev/null || true
                JOB_JSON=$(dx describe "$BASELINE_JOB" --json 2>/dev/null || echo "{}")
                GFA_FILE_ID=$(echo "$JOB_JSON" | python3 -c "
import sys,json
d=json.load(sys.stdin)
o=d.get('output',{})
f=o.get('gfa',{})
if isinstance(f,dict): print(f.get('\$dnanexus_link',''))
elif f: print(str(f))
" 2>/dev/null || echo "")
                [ -n "$GFA_FILE_ID" ] && dx download "$GFA_FILE_ID" -o "$BASELINE_GFA" 2>/dev/null || true
            fi
        fi
    fi
fi
# Fallback: local execution if DNAnexus unavailable
if [ ! -f "$BASELINE_GFA" ]; then
    echo "  Local baseline execution..."
    bash pipeline/baseline/build_baseline.sh "$INPUT_FASTA" "results/baseline"
fi
[ ! -f "$BASELINE_GFA" ] && { echo "FATAL: no baseline graph"; exit 1; }
echo "  Baseline: $BASELINE_GFA"

echo "[6/9] Chunks"
python3 pipeline/parallel/make_chunks.py
python3 pipeline/parallel/build_all_chunks.py

echo "[6b/9] Launching parallel chunk PGGB on DNAnexus"
bash dnanexus/run_parallel_chunks.sh

echo "[7/9] Merge (diagnostic disjoint union only — NOT_IMPLEMENTED)"
python3 pipeline/merge/merge_graphs.py
MERGED_GFA="results/merge/merged.gfa"

echo "[8/9] Benchmark"
python3 pipeline/benchmark/graph_stats.py
python3 pipeline/benchmark/compare_paths.py
bash pipeline/benchmark/benchmark_variants.sh 2>/dev/null || true
python3 pipeline/benchmark/build_report.py

echo "[9/9] Web JSON — compact bounded JSON only, no large GFA copy"
python3 pipeline/export/gfa_to_json.py "$BASELINE_GFA" --output "web/public/data/baseline.json" --label "baseline" 2>/dev/null || true
[ -f "$MERGED_GFA" ] && python3 pipeline/export/gfa_to_json.py "$MERGED_GFA" --output "web/public/data/merged_diagnostic.json" --label "merged_diagnostic" 2>/dev/null || true
python3 scripts/sync_web_results.py 2>/dev/null || true

echo ""
echo "=== Pipeline Summary ==="
echo "  baseline: PASS"
echo "  chunks: PASS"
echo "  parallel_execution: PASS"
echo "  stitch: NOT_IMPLEMENTED"
echo "  equivalence: NOT_RUN"
echo "  overall_status: PARTIAL"
echo "  (Overlap-aware stitching not yet implemented. Merged graph is diagnostic disjoint union only.)"
echo ""
echo "  Baseline: $BASELINE_GFA"
echo "  Merged (diagnostic): ${MERGED_GFA:-N/A}"

if [ "$UPLOAD" = "--upload" ]; then
    dx upload results/merge/merged.gfa --destination /graphs/merged/ 2>/dev/null || true
    dx upload results/baseline/baseline.gfa --destination /graphs/baseline/ 2>/dev/null || true
    dx upload results/benchmark/report.json --destination /benchmark/ 2>/dev/null || true
fi