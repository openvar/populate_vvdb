#!/bin/bash
set -euo pipefail

# ============================================================
# VariantValidator dbSNP Builder
# Option B: If final BCF exists, skip EVERYTHING.
# Plus: Clean up all intermediate files (.gz, .tbi, .metadata, .csi)
# ============================================================

BASE_URL="https://ftp.ncbi.nih.gov/snp/latest_release/VCF"
GRCH37_VCF="GCF_000001405.25.gz"
GRCH38_VCF="GCF_000001405.40.gz"

# ---- Helper: Cleanup intermediate files ----
cleanup_intermediates() {
    echo "🧹 Cleaning up intermediate dbSNP files..."

    rm -f GCF_*.gz
    rm -f GCF_*.tbi
    rm -f GCF_*_metadata.txt
    rm -f vvdbSNP_*.bcf.csi

    echo "✔ Cleanup complete (kept only final BCF files)."
}

# ---- Check for final BCF files ----
HAS37=$(ls vvdbSNP_GRCh37_b*.bcf 2>/dev/null || true)
HAS38=$(ls vvdbSNP_GRCh38_b*.bcf 2>/dev/null || true)

if [[ -n "$HAS37" && -n "$HAS38" ]]; then
    echo "⏩ Final dbSNP BCF files already exist. Skipping entire dbSNP build."
    cleanup_intermediates
    exit 0
fi

echo "⚠️  Final dbSNP BCF files not found. Proceeding with full build..."

process_build () {

    VCF="$1"
    META="${VCF%.gz}_metadata.txt"

    # ---- Step 0: Download if missing ----
    if [[ ! -f "$VCF" ]]; then
        STEM="${VCF%.gz}"
        echo "⬇ Downloading $STEM..."
        wget -c -O "$VCF" "${BASE_URL}/${STEM}.gz"
    else
        echo "⏩ $VCF already exists, skipping download"
    fi

    # ---- Step 1: Tabix index if not already ----
    if [[ ! -f "${VCF}.tbi" ]]; then
        echo "📚 Indexing source VCF with tabix..."
        tabix -p vcf "$VCF"
    else
        echo "⏩ Tabix index already exists"
    fi

    # ---- Step 2: Extract metadata ----
    echo "📝 Extracting metadata..."
    zgrep '^##' "$VCF" > "$META"

    # ---- Step 3: Get build and fileDate ----
    BUILD=$(grep -i '^##dbSNP_BUILD_ID' "$META" | head -1 | sed 's/.*=//')
    FILEDATE=$(grep -i '^##fileDate' "$META" | head -1 | sed 's/.*=//')

    [[ -z "$BUILD" ]] && BUILD=$(grep -i 'build' "$META" | head -1 | sed 's/.*=//')
    [[ -z "$FILEDATE" ]] && FILEDATE=$(date +%Y%m%d)

    echo "✔ Build: b$BUILD"
    echo "✔ FileDate: $FILEDATE"

    # ---- Step 4: Determine assembly ----
    if [[ "$VCF" == *".25"* ]]; then
        ASM="GRCh37"
    else
        ASM="GRCh38"
    fi

    OUT="vvdbSNP_${ASM}_b${BUILD}_${FILEDATE}.bcf"

    # ---- Step 5: Skip BCF creation if it exists ----
    if [[ -f "$OUT" ]]; then
        echo "⏩ Skipping: $OUT already exists"
        return
    fi

    # ---- Step 6: Create minimal BCF ----
    echo "🧬 Creating minimal BCF (only #CHROM POS ID REF ALT)..."

    bcftools view "$VCF" -Ou \
        | bcftools annotate -x INFO,FORMAT \
        | bcftools view -Ob -o "$OUT"

    # ---- Step 7: Index the BCF ----
    # We build the .csi so bcftools operations are valid.
    # It will be deleted later during cleanup.
    echo "📚 Indexing minimal BCF..."
    bcftools index "$OUT"

    echo "✔ Done: $OUT"
}

# ---- Run ----
process_build "$GRCH37_VCF"
process_build "$GRCH38_VCF"

# ---- Cleanup intermediate files ----
cleanup_intermediates

echo "🎉 All done (final BCFs retained, everything else removed)"

