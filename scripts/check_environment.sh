#!/usr/bin/env bash
# check_environment.sh
# Verify that required and optional tools are available.
# Exit code: 0 if all REQUIRED tools found, 1 otherwise.
# Prints a summary table.
#
# Tiers:
#   REQUIRED      — needed for the core genomics pipeline on DNAnexus
#   CONTAINER     — used inside Docker (pggb, vg, odgi — not on shell PATH)
#   WEB OPTIONAL  — needed only for the Next.js frontend (node/npm)
#   OPTIONAL      — nice-to-have but not blocking

set -uo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

all_required_found=true

check_tool() {
    local tool=$1
    local category=$2
    if command -v "$tool" &>/dev/null; then
        echo -e "${GREEN}FOUND${NC}    ${tool}"
        return 0
    else
        case "$category" in
            REQUIRED)
                echo -e "${RED}MISSING${NC}  ${tool} (REQUIRED — install to proceed)"
                all_required_found=false
                return 1
                ;;
            CONTAINER)
                echo -e "${YELLOW}IN-DOCKER${NC} ${tool} (available inside container, not on PATH)"
                return 0
                ;;
            WEB_OPTIONAL)
                echo -e "${CYAN}WEB-OPT${NC}  ${tool} (needed only for web UI development)"
                return 0
                ;;
            *)
                echo -e "${YELLOW}OPTIONAL${NC} ${tool} (not found)"
                return 0
                ;;
        esac
    fi
}

echo "=========================================="
echo "  Environment Check: Parallel Pangenome"
echo "=========================================="
echo ""

echo "--- Required Pipeline Tools ---"
check_tool "python3" "REQUIRED"
check_tool "pip3" "REQUIRED" || check_tool "pip" "REQUIRED"
check_tool "git" "REQUIRED"
check_tool "docker" "REQUIRED"

echo ""
echo "--- Cloud & Data Tools ---"
check_tool "dx" "REQUIRED"
check_tool "aws" "REQUIRED"

echo ""
echo "--- Bioinformatics Tools (installed natively) ---"
check_tool "samtools" "REQUIRED"
check_tool "bgzip" "REQUIRED"
check_tool "tabix" "REQUIRED"
check_tool "minimap2" "REQUIRED"
check_tool "bcftools" "REQUIRED"
check_tool "truvari" "REQUIRED"

echo ""
echo "--- Bioinformatics Tools (in Docker containers) ---"
check_tool "pggb" "CONTAINER"
check_tool "vg" "CONTAINER"
check_tool "odgi" "CONTAINER"

echo ""
echo "--- Web Development (optional for Vercel/frontend work) ---"
check_tool "node" "WEB_OPTIONAL"
check_tool "npm" "WEB_OPTIONAL"

echo ""
echo "--- Python Packages ---"
python3 -c "import yaml; print('FOUND    pyyaml')" 2>/dev/null \
    || echo -e "${YELLOW}MISSING${NC}  pyyaml (pip install pyyaml)"

echo ""
echo "=========================================="
if [ "$all_required_found" = true ]; then
    echo -e "${GREEN}All required pipeline tools are available.${NC}"
    exit 0
else
    echo -e "${RED}Some required pipeline tools are missing.${NC}"
    echo "Install them before proceeding with the HPRC pipeline."
    echo "See dnanexus/setup_workstation.sh for automated install."
    exit 1
fi
