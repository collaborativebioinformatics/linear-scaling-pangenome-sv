#!/usr/bin/env bash
# create_project_dirs.sh — Create DNAnexus project directory structure.
# Uses DX_PROJECT_CONTEXT_ID (automatically set inside DNAnexus Cloud Workstations)
# or falls back to DX_PROJECT_ID if set manually.
set -euo pipefail

# Auto-detect project context
PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"

echo "=== DNAnexus Project Directory Setup ==="
echo ""

if [ -z "$PROJECT_ID" ]; then
    echo "Not running inside DNAnexus (DX_PROJECT_CONTEXT_ID unset)."
    echo "Showing planned directory structure:"
    DRY_RUN=true
else
    echo "Using project: $PROJECT_ID"
    DRY_RUN=false
fi

DIRS=(
    /data/hprc
    /data/reference
    /data/prepared
    /graphs/baseline
    /graphs/chunks
    /graphs/merged
    /variants
    /benchmark
    /web
    /logs
)

for DIR in "${DIRS[@]}"; do
    if [ "$DRY_RUN" = true ]; then
        echo "  Would create: $DIR"
    else
        echo "  Creating $DIR ..."
        dx mkdir -p "$PROJECT_ID:$DIR" 2>/dev/null \
            && echo "    OK" \
            || echo "    FAILED (may already exist)"
    fi
done

echo ""
echo "=== Done ==="
if [ "$DRY_RUN" = true ]; then
    echo "Run inside a DNAnexus Cloud Workstation to actually create directories."
fi