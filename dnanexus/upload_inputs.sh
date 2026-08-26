#!/usr/bin/env bash
# upload_inputs.sh — Upload HPRC assemblies, reference, and config to DNAnexus.
# Run inside a DNAnexus Cloud Workstation after downloads complete.
set -euo pipefail

echo "=== Upload Inputs to DNAnexus ==="
echo ""

PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
if [ -z "$PROJECT_ID" ]; then
    echo "ERROR: Not running inside DNAnexus (DX_PROJECT_CONTEXT_ID unset)."
    exit 1
fi
echo "Project: $PROJECT_ID"
echo ""

# 1. HPRC assemblies
if ls work/downloads/*.fa 2>/dev/null; then
    echo "Uploading HPRC assemblies..."
    dx upload work/downloads/*.fa --destination /data/hprc/ --brief 2>/dev/null || true
    echo "  -> /data/hprc/"
else
    echo "  SKIP: no assemblies in work/downloads/"
fi
echo ""

# 2. Reference (if available)
if [ -f "reference/GRCh38.fa" ]; then
    echo "Uploading GRCh38 reference..."
    dx upload reference/GRCh38.fa --destination /data/reference/ --brief 2>/dev/null || true
elif [ -f "work/downloads/GRCh38.fa" ]; then
    echo "Uploading GRCh38 reference..."
    dx upload work/downloads/GRCh38.fa --destination /data/reference/ --brief 2>/dev/null || true
else
    echo "  SKIP: GRCh38.fa not found. Will fetch during preparation."
fi
echo ""

# 3. Config
echo "Uploading pipeline config..."
dx upload config/pipeline.yaml --destination /data/ --brief 2>/dev/null || true
dx upload config/samples.yaml --destination /data/ --brief 2>/dev/null || true
echo "  -> /data/"
echo ""

echo "=== Upload Complete ==="