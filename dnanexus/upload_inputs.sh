#!/usr/bin/env bash
# upload_inputs.sh — Upload prepared inputs to DNAnexus project.
# TODO (Ali): Add file-by-file upload with progress.
set -euo pipefail

echo "=== Upload Inputs to DNAnexus ==="
echo "REAL-DATA STEP: Requires dx login and selected project."
echo ""
echo "Run these commands inside DNAnexus Cloud Workstation:"
echo ""
echo "  # Upload HPRC assemblies"
echo "  dx upload work/downloads/*.fa --destination /data/hprc/"
echo ""
echo "  # Upload reference"
echo "  dx upload reference/GRCh38.fa --destination /data/reference/"
echo ""
echo "  # Upload config"
echo "  dx upload config/pipeline.yaml --destination /data/"