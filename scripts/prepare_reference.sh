#!/usr/bin/env bash
# prepare_reference.sh — Prepare GRCh38 chr21 reference.
#
# Staging order: local file → DNAnexus → external NCBI/Ensembl.
# Metadata calculated AFTER the final reference is established.
set -euo pipefail

echo "=== Prepare GRCh38 Reference (chr21) ==="

REF_DIR="work/reference"
META_DIR="results/preparation"
mkdir -p "$REF_DIR" "$META_DIR"

CHR21_FA="$REF_DIR/GRCh38_chr21.fa"
ACCESSION="GCF_000001405.40"
DNANEXUS_FILE_ID=""
SOURCE=""

# Step 1: Check local file
if [ -f "$CHR21_FA" ] && [ -s "$CHR21_FA" ]; then
    echo "  Reference already exists locally"
    SOURCE="local"
# Step 2: DNAnexus persistent storage
elif [ -n "${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}" ]; then
    PROJECT_ID="${DX_PROJECT_CONTEXT_ID:-${DX_PROJECT_ID:-}}"
    echo "  Checking DNAnexus /data/reference/..."
    DNA_FILE_ID=$(dx find data --name "GRCh38_chr21.fa" --path "$PROJECT_ID:/data/reference" --brief 2>/dev/null | head -1 || echo "")
    if [ -n "$DNA_FILE_ID" ]; then
        echo "  Staging from DNAnexus: $DNA_FILE_ID"
        dx download "$DNA_FILE_ID" -o "$CHR21_FA" 2>/dev/null && SOURCE="dnanexus" && DNANEXUS_FILE_ID="$DNA_FILE_ID"
    fi
fi

# Step 3: External NCBI/Ensembl fallback
if [ -z "$SOURCE" ]; then
    echo "  Downloading from external source..."
    NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_stuff/chr21.fa.gz"
    ENSEMBL_URL="https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
    if wget -q --timeout=30 "$NCBI_URL" -O "$REF_DIR/chr21.fa.gz" 2>/dev/null; then
        echo "  Downloaded from NCBI"; SOURCE="ncbi"
    elif wget -q --timeout=30 "$ENSEMBL_URL" -O "$REF_DIR/chr21.fa.gz" 2>/dev/null; then
        echo "  Downloaded from Ensembl"; SOURCE="ensembl"
    else
        echo "ERROR: Could not download GRCh38 chr21."; exit 1
    fi
    gunzip -f "$REF_DIR/chr21.fa.gz"; mv "$REF_DIR/chr21.fa" "$CHR21_FA"
    sed -i '' 's/^>.*/>GRCh38#0#chr21/' "$CHR21_FA" 2>/dev/null || sed -i 's/^>.*/>GRCh38#0#chr21/' "$CHR21_FA"
fi

# Calculate metadata AFTER reference is finalized
BP=$(grep -v '^>' "$CHR21_FA" | tr -d '\n' | wc -c)
SHA256=$(sha256sum "$CHR21_FA" 2>/dev/null | awk '{print $1}' || sha256 -r "$CHR21_FA" 2>/dev/null | awk '{print $1}' || echo "N/A")
HEADER=$(grep '^>' "$CHR21_FA" | head -1)

cat > "$META_DIR/reference_metadata.json" << JSONEOF
{
  "source": "$SOURCE",
  "dnanexus_file_id": "$DNANEXUS_FILE_ID",
  "assembly_accession": "$ACCESSION",
  "assembly_name": "GRCh38.p14",
  "chromosome": "chr21",
  "file": "$CHR21_FA",
  "length_bp": $BP,
  "sha256": "$SHA256",
  "normalized_header": "$HEADER"
}
JSONEOF
echo "Source: $SOURCE, Length: $BP bp, SHA256: $SHA256"
echo "=== Reference ready ==="
