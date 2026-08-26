#!/usr/bin/env bash
# check_environment.sh
# Verify that required and optional tools are available.
# Exit code: 0 if all REQUIRED tools found, 1 otherwise.
# Prints a summary table.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

check_tool() {
    local tool=$1
    local category=$2  # REQUIRED or OPTIONAL
    if command -v "$tool" &>/dev/null; then
        echo -e "${GREEN}FOUND${NC}    ${tool}"
        return 0
    else
        if [ "$category" = "REQUIRED" ]; then
            echo -e "${RED}MISSING${NC}  ${tool} (REQUIRED)"
            return 1
        else
            echo -e "${YELLOW}OPTIONAL${NC} ${tool} (not found)"
            return 0
        fi
    fi
}

echo "=========================================="
echo "  Environment Check: Parallel Pangenome"
echo "=========================================="
echo ""

all_required_found=true

# --- REQUIRED ---
echo "--- Required Tools ---"
check_tool "python3" "REQUIRED" || all_required_found=false
check_tool "git" "REQUIRED" || all_required_found=false
check_tool "node" "REQUIRED" || all_required_found=false
check_tool "npm" "REQUIRED" || all_required_found=false

echo ""
echo "--- Container Runtime ---"
check_tool "docker" "OPTIONAL"

echo ""
echo "--- Cloud Tools ---"
check_tool "dx" "OPTIONAL"
check_tool "aws" "OPTIONAL"

echo ""
echo "--- Bioinformatics Tools ---"
check_tool "samtools" "OPTIONAL"
check_tool "bgzip" "OPTIONAL"
check_tool "tabix" "OPTIONAL"
check_tool "minimap2" "OPTIONAL"
check_tool "pggb" "OPTIONAL"
check_tool "vg" "OPTIONAL"
check_tool "odgi" "OPTIONAL"
check_tool "bcftools" "OPTIONAL"
check_tool "truvari" "OPTIONAL"

echo ""
echo "--- Python Packages (core) ---"
python3 -c "import yaml; print('FOUND    pyyaml')" 2>/dev/null && true || echo -e "${YELLOW}OPTIONAL${NC} pyyaml (recommended)"

echo ""
echo "=========================================="
if [ "$all_required_found" = true ]; then
    echo -e "${GREEN}All required tools are available.${NC}"
    exit 0
else
    echo -e "${RED}Some required tools are missing.${NC}"
    echo "Install missing tools before proceeding."
    exit 1
fi
