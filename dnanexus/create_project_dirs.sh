#!/usr/bin/env bash
# create_project_dirs.sh — Create DNAnexus project directory structure.
# TODO (Ali): Add dx mkdir commands once project ID is configured.
set -euo pipefail

PROJECT_ID="${DX_PROJECT_ID:-}"

echo "DNAnexus Project Directory Setup"
echo "================================"

if [ -z "$PROJECT_ID" ]; then
    echo "DX_PROJECT_ID not set. Showing planned structure:"
fi

for DIR in \
    /data/hprc \
    /data/reference \
    /data/prepared \
    /graphs/baseline \
    /graphs/chunks \
    /graphs/merged \
    /variants \
    /benchmark \
    /web \
    /logs; do
    if [ -n "$PROJECT_ID" ]; then
        echo "  Creating $DIR ..."
        dx mkdir -p "$PROJECT_ID:$DIR" 2>/dev/null || echo "  FAILED: $DIR"
    else
        echo "  Would create: $DIR"
    fi
done

echo ""
echo "To enable: export DX_PROJECT_ID=project-xxxx"