#!/usr/bin/env bash
# setup_workstation.sh — Bootstrap a DNAnexus Cloud Workstation (Ubuntu 24.04)
# for the pangenome-parallel pipeline.
#
# Run ONCE inside the workstation after launching:
#   dx run --instance-type mem3_ssd1_v2_x16 --ssh app-cloud_workstation
#   git clone <repo> && cd pangenome-parallel
#   bash dnanexus/setup_workstation.sh
#
set -euo pipefail

echo "=== DNAnexus Workstation Setup ==="
echo "  OS: Ubuntu 24.04 (Cloud Workstation)"
python3 --version 2>/dev/null || echo "  python3 not found"
docker --version 2>/dev/null || echo "  docker not found"
dx --version 2>/dev/null || echo "  dx not found"
echo ""

# ── 1. System dependencies ──────────────────────────────────────────────────
echo "[1/4] Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y -qq \
    curl wget git \
    python3-pip python3-dev \
    samtools bcftools tabix \
    minimap2 build-essential \
    2>&1 | tail -1
echo "  done"

# ── 2. AWS CLI (for public HPRC S3 access) ──────────────────────────────────
echo "[2/4] Installing AWS CLI..."
if ! command -v aws &>/dev/null; then
    curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip"
    unzip -q /tmp/awscliv2.zip -d /tmp/aws
    sudo /tmp/aws/aws/install --update 2>&1 | tail -1
    rm -rf /tmp/aws /tmp/awscliv2.zip
fi
aws --version 2>&1 | head -1
echo "  done"

# ── 3. Python packages ──────────────────────────────────────────────────────
echo "[3/4] Installing Python packages..."
python3 -m pip install --quiet --upgrade pip wheel setuptools
python3 -m pip install --quiet pyyaml pytest awscli 2>&1 | tail -1
echo "  done"

# ── 4. Verify environment ───────────────────────────────────────────────────
# ── 5. Pull Docker images ───────────────────────────────────────────────────
echo "[5/5] Pulling Docker images..."
docker pull ghcr.io/pangenome/pggb:latest 2>&1 | tail -1 || true
docker pull quay.io/vgteam/vg:v1.74.1 2>&1 | tail -1 || true
echo "  done"
echo "[4/4] Verifying environment..."
echo ""
bash scripts/check_environment.sh || true
echo ""

# ── Done ────────────────────────────────────────────────────────────────────
echo "=== Setup Complete ==="
echo ""
echo "Next steps:"
echo "  1. python3 scripts/fetch_hprc_index.py"
echo "  2. python3 scripts/download_hprc.py --execute"
echo "  3. dx upload work/downloads/*.fa --destination /data/hprc/"
echo "  4. bash dnanexus/create_project_dirs.sh"
echo ""
echo "See dnanexus/README.md for the full pipeline."
