#!/bin/bash
set -euo pipefail

# ---- Settings ----
BASE_URL="https://ftp.ncbi.nih.gov/snp/latest_release/VCF"
GRCH37_VCF="GCF_000001405.25.gz"
GRCH38_VCF="GCF_000001405.40.gz"

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

    # ---- Step 5: Skip if already exists ----
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
    echo "📚 Indexing minimal BCF..."
    bcftools index "$OUT"

    echo "✔ Done: $OUT"
}

# ---- Run ----
process_build "$GRCH37_VCF"
process_build "$GRCH38_VCF"

echo "🎉 All done"