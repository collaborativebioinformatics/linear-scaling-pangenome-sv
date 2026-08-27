#!/usr/bin/env bash
# run_svim_asm.sh — Alternative assembly-based SV caller (SVIM-asm).
#
# Second opinion on dipcall. dipcall derives variants from paftools.js
# `call` over minimap2 alignments; SVIM-asm has its own signature-clustering
# logic and reports SV types dipcall does not (notably inversions and
# tandem/interspersed duplications). Disagreement between the two is a
# useful flag for representation-sensitive loci.
#
# Requires: svim-asm, minimap2, samtools.
#   conda install -c bioconda svim-asm minimap2 samtools
#
# Usage:
#   bash pipeline/linear/run_svim_asm.sh
#   bash pipeline/linear/run_svim_asm.sh HG00673
#
# Env overrides:
#   REF, OUTDIR, THREADS, DOWNLOAD_DIR, MANIFEST, MIN_SVLEN
set -euo pipefail

# --- config-derived defaults -----------------------------------------------
# The team's config/pipeline.yaml defines target.reference (e.g. GRCh38) and
# target.chromosome (e.g. chr21). scripts/prepare_reference.sh writes
# work/reference/<reference>_<chromosome>.fa. Deriving both here means these
# scripts keep working when the target region changes, instead of silently
# pointing at a chr21 file that no longer exists.
CONFIG="${CONFIG:-config/pipeline.yaml}"
_cfg() {  # _cfg <dotted.key> <fallback>
    if [ -f "$CONFIG" ] && command -v python3 &>/dev/null; then
        python3 -c "
import sys,yaml
try:
    c=yaml.safe_load(open('$CONFIG')) or {}
    for k in '$1'.split('.'): c=c[k]
    print(c)
except Exception:
    print('$2')
" 2>/dev/null || echo "$2"
    else
        echo "$2"
    fi
}
REF_NAME="${REF_NAME:-$(_cfg target.reference GRCh38)}"
CHROM="${CHROM:-$(_cfg target.chromosome chr21)}"

REF="${REF:-work/reference/${REF_NAME}_${CHROM}.fa}"
OUTDIR="${OUTDIR:-results/linear/svim_asm}"
DOWNLOAD_DIR="${DOWNLOAD_DIR:-work/downloads}"
MANIFEST="${MANIFEST:-work/manifests/hprc_selected.csv}"
MIN_SVLEN="${MIN_SVLEN:-50}"

if command -v nproc &>/dev/null; then
    THREADS="${THREADS:-$(( $(nproc) < 8 ? $(nproc) : 8 ))}"
else
    THREADS="${THREADS:-4}"
fi

echo "=== SVIM-asm assembly-based SV calling ==="

missing=0
for tool in svim-asm minimap2 samtools; do
    if ! command -v "$tool" &>/dev/null; then
        echo "  MISSING: $tool"
        missing=1
    fi
done
if [ "$missing" -ne 0 ]; then
    echo ""
    echo "  SKIP: SVIM-asm toolchain not installed."
    echo "  Install with: conda install -c bioconda svim-asm minimap2 samtools"
    exit 0
fi

if [ ! -f "$REF" ]; then
    echo "  FATAL: reference not found: $REF" >&2
    exit 1
fi
if [ ! -f "$MANIFEST" ]; then
    echo "  FATAL: manifest not found: $MANIFEST" >&2
    exit 1
fi

mkdir -p "$OUTDIR" work/linear/bam results/logs


# Resolve a manifest column BY HEADER NAME, not by position. The manifest
# column order is set by scripts/fetch_hprc_index.py; a positional awk here
# would break silently the first time anyone reorders it.
manifest_field() {  # manifest_field <sample_id> <haplotype_label> <column>
    python3 - "$MANIFEST" "$1" "$2" "$3" <<'PYMF'
import csv, sys
path, sample, label, col = sys.argv[1:5]
with open(path, newline="") as f:
    for r in csv.DictReader(f):
        if r.get("sample_id") == sample and r.get("haplotype_label") == label:
            print(r.get(col, ""))
            break
PYMF
}

find_assembly() {
    local name="$1"
    for ext in .fa.gz .fa .fasta .fna; do
        for dir in "$DOWNLOAD_DIR" /data/hprc; do
            if [ -f "${dir}/${name}${ext}" ]; then
                echo "${dir}/${name}${ext}"
                return 0
            fi
        done
    done
    return 1
}

# Align one haplotype assembly to the reference and produce a sorted BAM.
# Flags follow the SVIM-asm docs: asm5 for intra-species assembly-to-ref,
# --cs for the difference string, -r2k for a longer bandwidth so large
# indels are not split across records.
align_hap() {
    local fa="$1" tag="$2"
    local bam="work/linear/bam/${tag}.sorted.bam"
    if [ -f "${bam}.bai" ]; then
        echo "$bam"
        return 0
    fi
    minimap2 -a -x asm5 --cs -r2k -t "$THREADS" "$REF" "$fa" 2>/dev/null \
        | samtools sort -@ 2 -m 4G -o "$bam" -
    samtools index "$bam"
    echo "$bam"
}

if [ $# -ge 1 ]; then
    SAMPLES=("$@")
else
    # `mapfile` is bash 4+; macOS ships bash 3.2, so read the list portably.
    SAMPLES=()
    while IFS= read -r _s; do
        [ -n "$_s" ] && SAMPLES+=("$_s")
    done < <(tail -n +2 "$MANIFEST" | cut -d, -f1 | sort -u)
fi

if [ "${#SAMPLES[@]}" -eq 0 ]; then
    echo "  FATAL: no samples resolved from $MANIFEST" >&2
    exit 1
fi

echo "  Reference: $REF"
echo "  Threads:   $THREADS"
echo "  Samples:   ${SAMPLES[*]}"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    echo "--- $SAMPLE ---"

    MAT_NAME=$(manifest_field "$SAMPLE" maternal assembly_name)
    PAT_NAME=$(manifest_field "$SAMPLE" paternal assembly_name)

    if [ -z "$MAT_NAME" ] || [ -z "$PAT_NAME" ]; then
        echo "  SKIP: could not resolve both haplotypes for $SAMPLE" >&2
        continue
    fi
    if ! MAT_FA=$(find_assembly "$MAT_NAME"); then
        echo "  SKIP: assembly not downloaded: $MAT_NAME" >&2
        continue
    fi
    if ! PAT_FA=$(find_assembly "$PAT_NAME"); then
        echo "  SKIP: assembly not downloaded: $PAT_NAME" >&2
        continue
    fi

    echo "  Aligning hap1 (paternal)..."
    HAP1_BAM=$(align_hap "$PAT_FA" "${SAMPLE}_hap1")
    echo "  Aligning hap2 (maternal)..."
    HAP2_BAM=$(align_hap "$MAT_FA" "${SAMPLE}_hap2")

    SAMPLE_OUT="${OUTDIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUT"

    echo "  Running svim-asm diploid..."
    svim-asm diploid \
        --sample "$SAMPLE" \
        --min_sv_size "$MIN_SVLEN" \
        "$SAMPLE_OUT" "$HAP1_BAM" "$HAP2_BAM" "$REF" \
        2>&1 | tee "results/logs/svim_asm_${SAMPLE}.log"

    RAW_VCF="${SAMPLE_OUT}/variants.vcf"
    if [ ! -f "$RAW_VCF" ]; then
        echo "  ERROR: expected $RAW_VCF was not produced" >&2
        exit 1
    fi

    if command -v bcftools &>/dev/null; then
        SV_VCF="${OUTDIR}/${SAMPLE}.sv.vcf.gz"
        bcftools sort -Oz -o "$SV_VCF" "$RAW_VCF"
        bcftools index -t "$SV_VCF"
        N_SV=$(bcftools index -n "$SV_VCF")
        echo "  SV calls: $SV_VCF  ($N_SV records)"
    else
        echo "  Raw calls: $RAW_VCF"
    fi
    echo ""
done

echo "=== SVIM-asm complete -> $OUTDIR ==="
echo ""
echo "NOTE: SVIM-asm emits symbolic ALT alleles (<DEL>, <INV>, <DUP:TANDEM>)."
echo "For Truvari comparison against sequence-resolved graph calls, pass"
echo "  --reference \$REF --dup-to-ins"
echo "so symbolic records are resolved and DUPs are matched as insertions."
