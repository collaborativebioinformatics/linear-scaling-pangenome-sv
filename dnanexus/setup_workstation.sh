#!/usr/bin/env bash
# setup_workstation.sh -- Bootstrap a DNAnexus Cloud Workstation (Ubuntu 24.04)
# for the pangenome-parallel pipeline.
set -euo pipefail

echo "=== DNAnexus Workstation Setup ==="
python3 --version 2>/dev/null || echo "  python3 not found"
docker --version 2>/dev/null || echo "  docker not found"
dx --version 2>/dev/null || echo "  dx not found"
echo ""

export PATH="$HOME/.local/bin:$PATH"
if ! echo ":$PATH:" | grep -q ":$HOME/.local/bin:"; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
fi

echo "[1/5] Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y -qq \
    curl wget git unzip \
    python3-pip python3-dev python3-venv \
    samtools bcftools tabix \
    minimap2 build-essential \
    2>&1 | tail -1
echo "  done"

echo "[2/5] Installing AWS CLI..."
if ! command -v aws &>/dev/null; then
    curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip"
    unzip -q /tmp/awscliv2.zip -d /tmp/aws
    sudo /tmp/aws/aws/install --update 2>&1 | tail -1
    rm -rf /tmp/aws /tmp/awscliv2.zip
fi
aws --version 2>&1 | head -1
echo "  done"

echo "[3/5] Setting up Python environment..."
if [ ! -d "$HOME/.pipeline-env" ]; then
    python3 -m venv "$HOME/.pipeline-env"
fi
source "$HOME/.pipeline-env/bin/activate"
pip install --quiet --upgrade pip wheel setuptools 2>&1 | tail -1
pip install --quiet pyyaml pytest awscli 2>&1 | tail -1
echo "  Virtual env: $HOME/.pipeline-env"

echo "  Installing Truvari..."
pip install --quiet Truvari 2>&1 | tail -1 || true
if ! command -v truvari &>/dev/null; then
    python3 -m pip install --user Truvari 2>&1 | tail -1 || true
fi
truvari --version 2>/dev/null || echo "  Truvari: check manually"
echo "  done"

export PATH="$HOME/.pipeline-env/bin:$PATH"
echo 'source "$HOME/.pipeline-env/bin/activate" 2>/dev/null || true' >> ~/.bashrc

echo "[4/5] Pulling Docker images..."
docker pull ghcr.io/pangenome/pggb:latest 2>&1 | tail -1 || true
docker pull quay.io/vgteam/vg:v1.74.1 2>&1 | tail -1 || true
echo "  done"

echo "[5/5] Verifying environment..."
bash scripts/check_environment.sh || true
echo ""
echo "=== Setup Complete ==="
echo ""
echo "Next steps:"
echo "  source ~/.bashrc"
echo "  1. python3 scripts/fetch_hprc_index.py"
echo "  2. python3 scripts/download_hprc.py --execute"
echo "  3. dx upload work/downloads/*.fa.gz --destination /data/hprc/"
echo "  4. bash dnanexus/create_project_dirs.sh"
