#!/usr/bin/env bash
# stage_inputs.sh — Stage HPRC assemblies from DNAnexus project storage.
#
# Inspects Group11_2026:/data/hprc/ using dx commands.
# Downloads existing files (.fa.gz + .fai + .gzi + .md5) to work/downloads/.
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

# Required assemblies and their sidecar files
declare -a ASSEMBLY_NAMES
ASSEMBLY_NAMES=(
    "HG00673_mat_hprc_r2_v1.0.1"
    "HG00673_pat_hprc_r2_v1.0.1"
    "HG00733_mat_hprc_r2_v1.0.1"
    "HG00733_pat_hprc_r2_v1.0.1"
)

HPRC_REMOTE_DIR="/data/hprc"
LOCAL_DIR="work/downloads"
mkdir -p "$LOCAL_DIR"

# Track which assemblies are missing from DNAnexus
declare -a MISSING_FROM_DNANEXUS

for NAME in "${ASSEMBLY_NAMES[@]}"; do
    BASE_FILE="${NAME}.fa.gz"
    echo -n "  $BASE_FILE ... "

    # Check if already staged locally
    if [ -f "$LOCAL_DIR/$BASE_FILE" ] && [ -s "$LOCAL_DIR/$BASE_FILE" ]; then
        SIZE_MB=$(du -h "$LOCAL_DIR/$BASE_FILE" 2>/dev/null | cut -f1)
        echo "FOUND IN DNANEXUS (already staged, ${SIZE_MB:-?})"
        continue
    fi

    # Check DNAnexus project storage using --name flag (correct syntax)
    DNANEXUS_RESULT=$(dx find data --name "$BASE_FILE" --path "$PROJECT_ID:$HPRC_REMOTE_DIR" --brief 2>/dev/null | head -1 || echo "")

    if [ -n "$DNANEXUS_RESULT" ]; then
        echo "FOUND IN DNANEXUS -> staging..."
        # Download the main .fa.gz file
        dx download "$PROJECT_ID:$HPRC_REMOTE_DIR/$BASE_FILE" -o "$LOCAL_DIR/$BASE_FILE" 2>/dev/null
        if [ -f "$LOCAL_DIR/$BASE_FILE" ] && [ -s "$LOCAL_DIR/$BASE_FILE" ]; then
            SIZE_MB=$(du -h "$LOCAL_DIR/$BASE_FILE" 2>/dev/null | cut -f1)
            echo "    STAGED: $BASE_FILE ($SIZE_MB)"

            # Stage sidecar files if they exist on DNAnexus
            for SIDECAR in ".fai" ".gzi" ".md5"; do
                SIDECAR_FILE="${NAME}.fa.gz${SIDECAR}"
                SIDECAR_EXISTS=$(dx find data --name "$SIDECAR_FILE" --path "$PROJECT_ID:$HPRC_REMOTE_DIR" --brief 2>/dev/null | head -1 || echo "")
                if [ -n "$SIDECAR_EXISTS" ]; then
                    dx download "$PROJECT_ID:$HPRC_REMOTE_DIR/$SIDECAR_FILE" -o "$LOCAL_DIR/$SIDECAR_FILE" 2>/dev/null
                    echo "    STAGED: $SIDECAR_FILE"
                fi
            done

            # Validate gzip integrity if gzip is available
            if command -v gzip &>/dev/null; then
                if gzip -t "$LOCAL_DIR/$BASE_FILE" 2>/dev/null; then
                    echo "    INTEGRITY OK"
                else
                    echo "    WARNING: gzip integrity check FAILED"
                fi
            fi

            # Check MD5 if .md5 file was staged
            MD5_FILE="$LOCAL_DIR/${NAME}.fa.gz.md5"
            if [ -f "$MD5_FILE" ]; then
                EXPECTED=$(cat "$MD5_FILE" | tr -d ' \t\n')
                ACTUAL=$(md5sum "$LOCAL_DIR/$BASE_FILE" 2>/dev/null | cut -d' ' -f1 || \
                         md5 -r "$LOCAL_DIR/$BASE_FILE" 2>/dev/null | cut -d' ' -f1 || echo "")
                if [ -n "$ACTUAL" ] && [ "$EXPECTED" = "$ACTUAL" ]; then
                    echo "    MD5 OK"
                else
                    echo "    WARNING: MD5 mismatch (expected=$EXPECTED, actual=$ACTUAL)"
                fi
            fi
        else
            echo "    STAGED FAILED"
            MISSING_FROM_DNANEXUS+=("$BASE_FILE")
        fi
    else
        echo "MISSING"
        MISSING_FROM_DNANEXUS+=("$BASE_FILE")
    fi
done

echo ""

# Summary
echo "=== Stage Summary ==="
STAGED_COUNT=0
for NAME in "${ASSEMBLY_NAMES[@]}"; do
    if [ -f "$LOCAL_DIR/${NAME}.fa.gz" ] && [ -s "$LOCAL_DIR/${NAME}.fa.gz" ]; then
        STAGED_COUNT=$((STAGED_COUNT + 1))
    fi
done
echo "  Total required: ${#ASSEMBLY_NAMES[@]}"
echo "  Staged locally: $STAGED_COUNT"
echo "  Missing from DNAnexus: ${#MISSING_FROM_DNANEXUS[@]}"

if [ ${#MISSING_FROM_DNANEXUS[@]} -gt 0 ]; then
    echo ""
    echo "MISSING assemblies (will be downloaded from HPRC public S3):"
    for F in "${MISSING_FROM_DNANEXUS[@]}"; do
        echo "  - $F"
    done
    echo ""
    echo "Run: python3 scripts/download_hprc.py --execute"
    echo "Then: bash dnanexus/upload_inputs.sh"
    echo "Then re-run: bash dnanexus/stage_inputs.sh"
fi

echo ""
echo "=== Stage Complete ==="