#!/usr/bin/env bash
# run_parallel_chunks.sh — Run chunk graphs in parallel on DNAnexus.
# TODO (Ali/Quang): Convert to scatter-gather dx workflow with status monitoring.
set -euo pipefail

echo "=== Run Parallel Chunks on DNAnexus ==="
echo "REAL-DATA STEP: Requires project with uploaded data."
echo ""
echo "For each chunk in chunk_manifest.tsv:"
echo ""
echo '  dx run pggb \'
echo '    -i fasta=/data/prepared/chunk_0001.fa \'
echo '    --instance-type mem3_ssd1_v2_x16 \'
echo '    --name "Chunk 0001" \'
echo '    --destination /graphs/chunks/'
echo ""
echo "Then collect and merge:"
echo "  python3 pipeline/merge/merge_graphs.py"
echo "  dx upload results/merge/merged.gfa --destination /graphs/merged/"