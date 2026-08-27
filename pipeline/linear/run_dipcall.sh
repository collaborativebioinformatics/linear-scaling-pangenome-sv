#!/usr/bin/env bash
# run_dipcall.sh — Assembly-based SV calling with dipcall.
#
# Orthogonal validation for the graph pipeline: aligns each sample's two
# haplotype assemblies to GRCh38 with minimap2 (via dipcall's generated
# makefile) and emits a diploid VCF. Variants found here are independent
# of PGGB, so agreement is real evidence that graph-derived calls are not
# construction artifacts.
#
# Requires: dipcall (run-dipcall on PATH), minimap2, k8, samtools, htsbox.
#   conda install -c bioconda dipcall minimap2 k8 samtools
#
# Usage:
#   bash pipeline/linear/run_dipcall.sh                 # all samples in manifest
#   bash pipeline/linear/run_dipcall.sh HG00673         # one sample
#   REF=work/reference/GRCh38_chr21.fa THREADS=16 \
#     bash pipeline/linear/run_dipcall.sh HG00673
#
# Env overrides:
#   REF        reference FASTA           (default work/reference/GRCh38_chr21.fa)
#   OUTDIR     output directory          (default results/linear/dipcall)
#   THREADS    minimap2 threads          (default: nproc, capped at 8)
#   MIN_SVLEN  min |SVLEN| for SV subset (default 50)
#   MALE       set to 1 for male samples (adds -x hs38.PAR.bed handling)
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
OUTDIR="${OUTDIR:-results/linear/dipcall}"
DOWNLOAD_DIR="${DOWNLOAD_DIR:-work/downloads}"
MANIFEST="${MANIFEST:-work/manifests/hprc_selected.csv}"
MIN_SVLEN="${MIN_SVLEN:-50}"
MALE="${MALE:-0}"

if command -v nproc &>/dev/null; then
    THREADS="${THREADS:-$(( $(nproc) < 8 ? $(nproc) : 8 ))}"
else
    THREADS="${THREADS:-4}"
fi

echo "=== dipcall assembly-based SV calling ==="

# --- Preflight -------------------------------------------------------------
missing=0
for tool in run-dipcall minimap2 k8 samtools; do
    if ! command -v "$tool" &>/dev/null; then
        echo "  MISSING: $tool"
        missing=1
    fi
done
if [ "$missing" -ne 0 ]; then
    echo ""
    echo "  SKIP: dipcall toolchain not installed."
    echo "  Install with:"
    echo "    conda install -c bioconda dipcall minimap2 k8 samtools htslib"
    echo ""
    echo "  This is an OPTIONAL validation step; the graph pipeline does not"
    echo "  depend on it. Exiting 0 so the pipeline continues."
    exit 0
fi

if [ ! -f "$REF" ]; then
    echo "  FATAL: reference not found: $REF" >&2
    echo "  Run: bash scripts/prepare_reference.sh" >&2
    exit 1
fi
if [ ! -f "$MANIFEST" ]; then
    echo "  FATAL: manifest not found: $MANIFEST" >&2
    echo "  Run: python3 scripts/fetch_hprc_index.py" >&2
    exit 1
fi

mkdir -p "$OUTDIR" results/logs

# --- Resolve the assembly file for a given assembly_name -------------------

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

# --- Which samples? --------------------------------------------------------
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

    # Pull the mat/pat assembly_names for this sample from the manifest.
    MAT_NAME=$(manifest_field "$SAMPLE" maternal assembly_name)
    PAT_NAME=$(manifest_field "$SAMPLE" paternal assembly_name)

    if [ -z "$MAT_NAME" ] || [ -z "$PAT_NAME" ]; then
        echo "  SKIP: could not resolve both haplotypes for $SAMPLE" >&2
        continue
    fi

    if ! MAT_FA=$(find_assembly "$MAT_NAME"); then
        echo "  SKIP: assembly not downloaded: $MAT_NAME" >&2
        echo "        Run: python3 scripts/download_hprc.py --execute" >&2
        continue
    fi
    if ! PAT_FA=$(find_assembly "$PAT_NAME"); then
        echo "  SKIP: assembly not downloaded: $PAT_NAME" >&2
        continue
    fi

    PREFIX="${OUTDIR}/${SAMPLE}"
    MAKEFILE="${PREFIX}.mak"

    # dipcall convention: hap1 = paternal, hap2 = maternal.
    DIPCALL_ARGS=(-t "$THREADS")
    if [ "$MALE" = "1" ]; then
        # For male samples chrX/chrY need PAR handling; dipcall ships the BED.
        PAR_BED="${PAR_BED:-}"
        if [ -n "$PAR_BED" ] && [ -f "$PAR_BED" ]; then
            DIPCALL_ARGS+=(-x "$PAR_BED")
            echo "  Male sample: using PAR bed $PAR_BED"
        else
            echo "  WARNING: MALE=1 but PAR_BED unset; chrX/chrY calls may be wrong" >&2
        fi
    fi

    echo "  hap1 (pat): $PAT_FA"
    echo "  hap2 (mat): $MAT_FA"
    echo "  Generating makefile -> $MAKEFILE"

    run-dipcall "${DIPCALL_ARGS[@]}" "$PREFIX" "$REF" "$PAT_FA" "$MAT_FA" \
        > "$MAKEFILE"

    echo "  Running dipcall (make -j2)..."
    make -j2 -f "$MAKEFILE" 2>&1 | tee "results/logs/dipcall_${SAMPLE}.log"

    DIP_VCF="${PREFIX}.dip.vcf.gz"
    DIP_BED="${PREFIX}.dip.bed"

    if [ ! -f "$DIP_VCF" ]; then
        echo "  ERROR: expected $DIP_VCF was not produced" >&2
        exit 1
    fi

    echo "  Confident regions: $DIP_BED"
    echo "  Raw calls:         $DIP_VCF"

    # --- SV subset (>= MIN_SVLEN) for Truvari comparison -------------------
    if command -v bcftools &>/dev/null; then
        SV_VCF="${PREFIX}.sv.vcf.gz"
        # dipcall emits explicit REF/ALT (no symbolic alleles), so SV length
        # is |len(ALT) - len(REF)|. Keep biallelic records only.
        bcftools view -f PASS "$DIP_VCF" \
          | bcftools norm -m-any -f "$REF" \
          | bcftools filter -i "abs(strlen(ALT)-strlen(REF)) >= ${MIN_SVLEN}" \
          | bcftools sort -Oz -o "$SV_VCF"
        bcftools index -t "$SV_VCF"
        N_SV=$(bcftools index -n "$SV_VCF")
        echo "  SV subset (>=${MIN_SVLEN}bp): $SV_VCF  ($N_SV records)"
    else
        echo "  NOTE: bcftools not found; skipping SV subset." >&2
    fi
    echo ""
done

echo "=== dipcall complete -> $OUTDIR ==="
echo "Next: bash pipeline/benchmark/benchmark_variants.sh"
