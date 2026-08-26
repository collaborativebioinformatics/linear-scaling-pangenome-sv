#!/usr/bin/env bash
# download_results.sh — Download results from DNAnexus to local.
# TODO (Ali): Add selective download by file type/pattern.
set -euo pipefail

echo "=== Download Results from DNAnexus ==="
echo ""
echo "Run these commands locally:"
echo ""
echo "  # Download benchmark report"
echo "  dx download /benchmark/report.json"
echo ""
echo "  # Download merged graph"
echo "  dx download /graphs/merged/merged.gfa"
echo ""
echo "  # Download web data"
echo "  dx download /web/latest.json -o web/public/data/latest.json"
echo ""
echo "  # Sync to web"
echo "  python3 scripts/sync_web_results.py"