#!/usr/bin/env bash
# upload_inputs.sh — Upload HPRC assemblies (.fa.gz), reference, config and manifest to DNAnexus.
# Run inside a DNAnexus Cloud Workstation after downloads complete.
# This makes DNAnexus project storage the source of truth for input data.
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

# 1. HPRC assemblies (.fa.gz — the native HPRC format)
if ls work/downloads/*.fa.gz 2>/dev/null; then
    echo "Uploading HPRC assemblies (.fa.gz)..."
    for GZ in work/downloads/*.fa.gz; do
        NAME=$(basename "$GZ")
        # Check if already exists in DNAnexus storage
        EXISTING=$(dx find data --path "$PROJECT_ID:/data/hprc/$NAME" --brief 2>/dev/null | head -1 || echo "")
        if [ -n "$EXISTING" ]; then
            echo "  EXISTS: /data/hprc/$NAME (skipping upload)"
        else
            dx upload "$GZ" --destination /data/hprc/ --brief 2>/dev/null || true
            echo "  UPLOADED: /data/hprc/$NAME"
        fi
    done
    echo "  -> /data/hprc/"
else
    echo "  SKIP: no .fa.gz files in work/downloads/"
fi
echo ""

# 2. HPRC manifest
if [ -f "work/manifests/hprc_selected.csv" ]; then
    echo "Uploading HPRC manifest..."
    dx upload work/manifests/hprc_selected.csv --destination /data/hprc/ --brief 2>/dev/null || true
    echo "  -> /data/hprc/hprc_selected.csv"
else
    echo "  SKIP: manifest not found"
fi
echo ""

# 3. Reference (decompressed .fa for local use)
if [ -f "work/reference/GRCh38_chr21.fa" ]; then
    echo "Uploading GRCh38 chr21 reference..."
    EXISTING=$(dx find data --path "$PROJECT_ID:/data/reference/GRCh38_chr21.fa" --brief 2>/dev/null | head -1 || echo "")
    if [ -n "$EXISTING" ]; then
        echo "  EXISTS: /data/reference/GRCh38_chr21.fa (skipping upload)"
    else
        dx upload work/reference/GRCh38_chr21.fa --destination /data/reference/ --brief 2>/dev/null || true
        echo "  UPLOADED: /data/reference/"
    fi
else
    echo "  SKIP: GRCh38_chr21.fa not found"
fi
echo ""

# 4. Config
echo "Uploading pipeline config..."
dx upload config/pipeline.yaml --destination /data/ --brief 2>/dev/null || true
dx upload config/samples.yaml --destination /data/ --brief 2>/dev/null || true
echo "  -> /data/"
echo ""

echo "=== Upload Complete ==="