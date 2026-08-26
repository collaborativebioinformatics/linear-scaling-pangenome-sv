#!/usr/bin/env bash
# stage_inputs.sh — Stage HPRC assemblies from DNAnexus project storage.
#
# Inspects Group11_2026:/data/hprc/ using dx commands.
# Downloads existing files to work/downloads/.
# Reports FOUND IN DNANEXUS / STAGED / MISSING for every required assembly.
# Only falls back to HPRC public S3 for genuinely missing files.
#
# Usage:
#   bash dnanexus/stage_inputs.sh

set -euo pipefail

echo "=== Stage Inputs from DNAnexus Project Storage ==="
echo ""

PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
if [ -z "$PROJECT_ID" ]; then
    echo "ERROR: Not inside DNAnexus (DX_PROJECT_CONTEXT_ID unset)."
    exit 1
fi
echo "Project: $PROJECT_ID"
echo ""

# Required assemblies
REQUIRED_FILES=(
    "HG00673_mat_hprc_r2_v1.0.1.fa.gz"
    "HG00673_pat_hprc_r2_v1.0.1.fa.gz"
    "HG00733_mat_hprc_r2_v1.0.1.fa.gz"
    "HG00733_pat_hprc_r2_v1.0.1.fa.gz"
)

HPRC_REMOTE_DIR="/data/hprc"
LOCAL_DIR="work/downloads"

mkdir -p "$LOCAL_DIR"

# Track which files are missing from DNAnexus
declare -a MISSING_FROM_DNANEXUS

for FILE in "${REQUIRED_FILES[@]}"; do
    echo -n "  $FILE ... "

    # Check if already staged locally
    if [ -f "$LOCAL_DIR/$FILE" ] && [ -s "$LOCAL_DIR/$FILE" ]; then
        SIZE_MB=$(du -h "$LOCAL_DIR/$FILE" 2>/dev/null | cut -f1)
        echo "FOUND IN DNANEXUS (already staged, ${SIZE_MB:-?})"
        continue
    fi

    # Check DNAnexus project storage using dx
    DNANEXUS_RESULT=$(dx find data --path "$PROJECT_ID:$HPRC_REMOTE_DIR/$FILE" --brief 2>/dev/null | head -1 || echo "")

    if [ -n "$DNANEXUS_RESULT" ]; then
        # File exists in DNAnexus — download it
        echo "FOUND IN DNANEXUS -> staging..."
        dx download "$PROJECT_ID:$HPRC_REMOTE_DIR/$FILE" -o "$LOCAL_DIR/$FILE" 2>/dev/null
        if [ -f "$LOCAL_DIR/$FILE" ] && [ -s "$LOCAL_DIR/$FILE" ]; then
            SIZE_MB=$(du -h "$LOCAL_DIR/$FILE" 2>/dev/null | cut -f1)
            echo "    STAGED ($SIZE_MB)"
        else
            echo "    STAGED FAILED"
            MISSING_FROM_DNANEXUS+=("$FILE")
        fi
    else
        echo "MISSING"
        MISSING_FROM_DNANEXUS+=("$FILE")
    fi
done

echo ""

# Summary
echo "=== Stage Summary ==="
echo "  Total required: ${#REQUIRED_FILES[@]}"
STAGED_COUNT=0
for FILE in "${REQUIRED_FILES[@]}"; do
    if [ -f "$LOCAL_DIR/$FILE" ] && [ -s "$LOCAL_DIR/$FILE" ]; then
        STAGED_COUNT=$((STAGED_COUNT + 1))
    fi
done
echo "  Staged locally: $STAGED_COUNT"
echo "  Missing from DNAnexus: ${#MISSING_FROM_DNANEXUS[@]}"

if [ ${#MISSING_FROM_DNANEXUS[@]} -gt 0 ]; then
    echo ""
    echo "MISSING assemblies (will be downloaded from HPRC public S3):"
    for FILE in "${MISSING_FROM_DNANEXUS[@]}"; do
        echo "  - $FILE"
    done
    echo ""
    echo "Run: python3 scripts/download_hprc.py --execute"
    echo "Then: bash dnanexus/upload_inputs.sh"
    echo "Then re-run: bash dnanexus/stage_inputs.sh"
fi

echo ""
echo "=== Stage Complete ==="