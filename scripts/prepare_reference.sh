#!/usr/bin/env bash
# prepare_reference.sh — Download GRCh38 reference (chr21 only for smoke test).
#
# We download only chr21 to minimize storage. The full GRCh38 reference
# is available from multiple sources; we use the NCBI GRCh38 primary assembly.
#
# Output: work/reference/GRCh38_chr21.fa
#         results/preparation/reference_metadata.json

set -euo pipefail

echo "=== Prepare GRCh38 Reference (chr21) ==="

REF_DIR="work/reference"
META_DIR="results/preparation"
mkdir -p "$REF_DIR" "$META_DIR"

# Try multiple sources for chr21
CHR21_FA="$REF_DIR/GRCh38_chr21.fa"

if [ -f "$CHR21_FA" ] && [ -s "$CHR21_FA" ]; then
    SIZE_MB=$(du -h "$CHR21_FA" | cut -f1)
    echo "  Reference already exists: $CHR21_FA ($SIZE_MB)"
else
    echo "  Downloading GRCh38 chr21..."

    # Source 1: NCBI (preferred)
    NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_stuff/chr21.fa.gz"
    # Source 2: ENSEMBL fallback
    ENSEMBL_URL="https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"

    # Try NCBI first
    if wget -q --timeout=30 "$NCBI_URL" -O "$REF_DIR/chr21.fa.gz" 2>/dev/null; then
        echo "  Downloaded from NCBI"
    elif wget -q --timeout=30 "$ENSEMBL_URL" -O "$REF_DIR/chr21.fa.gz" 2>/dev/null; then
        echo "  Downloaded from Ensembl"
    else
        echo "  ERROR: Could not download GRCh38 chr21 from any source."
        echo "  Manual download:"
        echo "    wget $NCBI_URL -O $REF_DIR/chr21.fa.gz"
        echo "    gunzip $REF_DIR/chr21.fa.gz"
        echo "    mv $REF_DIR/chr21.fa $CHR21_FA"
        exit 1
    fi

    gunzip -f "$REF_DIR/chr21.fa.gz"
    mv "$REF_DIR/chr21.fa" "$CHR21_FA"

    # Normalize header
    sed -i '' 's/^>.*/>GRCh38#0#chr21/' "$CHR21_FA" 2>/dev/null || \
    sed -i 's/^>.*/>GRCh38#0#chr21/' "$CHR21_FA"

    echo "  Saved: $CHR21_FA"
fi

# Count bases
BP=$(grep -v '^>' "$CHR21_FA" | tr -d '\n' | wc -c)
echo "  Length: $BP bp"

# Write metadata
cat > "$META_DIR/reference_metadata.json" << JSONEOF
{
  "source": "GRCh38.p14",
  "chromosome": "chr21",
  "file": "$CHR21_FA",
  "length_bp": $BP,
  "normalized_header": "GRCh38#0#chr21"
}
JSONEOF
echo "  Metadata: $META_DIR/reference_metadata.json"
echo "  Metadata: $META_DIR/reference_metadata.json"

# Stage from DNAnexus project storage if available
STAGED=false
PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
if [ -n "$PROJECT_ID" ]; then
    DNA_FILE=$(dx find data --name "GRCh38_chr21.fa" --path "$PROJECT_ID:/data/reference" --brief 2>/dev/null | head -1 || echo "")
    if [ -n "$DNA_FILE" ]; then
        echo "  Staging from DNAnexus: $DNA_FILE"
        dx download "$PROJECT_ID:/data/reference/GRCh38_chr21.fa" -o "$CHR21_FA" 2>/dev/null && STAGED=true
    fi
fi
if [ "$STAGED" = true ]; then
    echo "  Staged from DNAnexus: $CHR21_FA"
fi

echo "=== Reference ready ==="
echo "=== Reference ready ==="