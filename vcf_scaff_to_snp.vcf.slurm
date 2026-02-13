#!/bin/bash
#SBATCH --job-name=gvcf2snp

set -euo pipefail

start=$(date +%s)
echo "Host: $HOSTNAME"

# --------------------
# Conda environment
# --------------------
dc=$(awk '{print $1}' "$2/conda.txt")
ec=$(awk '{print $2}' "$2/conda.txt")

source "${dc}/conda.sh"
conda activate "${ec}"

# --------------------
# Parameters
# --------------------
nind=$(wc -l < "$3")
minMQ=40
minDP=5

maxDP=$(awk '{print $2}' "$5"/*/*depth.txt | sort -nr | head -1 | awk '{print $1*3}')

minmaf=$(awk -v n="$nind" 'BEGIN {printf "%.6f\n", 2/(2*n)}')
minmafi=$(awk -v m="$minmaf" 'BEGIN {printf "%.6f\n", 1-m}')

sc1=$(awk '{print $1}' "$6")
sc2=$(awk '{print $2}' "$6")

TMP=${TMPDIR:-/tmp}

echo "Individuals: $nind"
echo "minMAF: $minmaf"
echo "maxDP: $maxDP"

# =====================================================
# 1. CONCATENATE SCAFFOLD VCFS
# =====================================================

if [[ -f "$1/gatk.all_indivs.vcf.gz" ]]; then
    echo "Concatenated VCF exists. Skipping."
else
    echo "Concatenating scaffold VCFs..."

    printf "%s\n" "$1"/scaff_vcfs/*.gatk.vcf.gz | sort -V > "$1/scaff_vcfs/vcf_list.txt"

    bcftools concat \
        -f "$1/scaff_vcfs/vcf_list.txt" \
        -Oz -o "$1/gatk.all_indivs.vcf.gz"
fi

bcftools index -f "$1/gatk.all_indivs.vcf.gz"
gatk IndexFeatureFile -I "$1/gatk.all_indivs.vcf.gz"

# Optional cleanup
if [[ "${7:-F}" == "T" ]]; then
    rm -f "$1"/indiv_scaff_gvcfs/*.rb.g.vcf.gz
    rm -f "$1"/scaff_vcfs/*.gatk.vcf.gz
fi

# =====================================================
# 2. SELECT VARIANT TYPES
# =====================================================

if [[ ! -f "$1/gatk.snp.indel.mnv.vcf.gz" ]]; then

    echo "Selecting SNP/INDEL/MNV..."

    gatk SelectVariants \
        -V "$1/gatk.all_indivs.vcf.gz" \
        -select-type SNP \
        -select-type INDEL \
        -select-type MNV \
        -O "$1/gatk.snp.indel.mnv.vcf.gz" \
        --tmp-dir "$TMP"
fi

bcftools index -f "$1/gatk.snp.indel.mnv.vcf.gz"

# =====================================================
# 3. SOFT FILTERING
# =====================================================

SOFT_TMP="$1/gatk.snp.indel.mnv.filtered_soft.tmp.vcf.gz"
SOFT_FINAL="$1/gatk.snp.indel.mnv.filtered_soft.vcf.gz"
HARD_FINAL="$1/gatk.snp.indel.mnv.filtered_hard.vcf.gz"

if [[ ! -f "$SOFT_TMP" && ! -f "$SOFT_FINAL" ]]; then

    echo "Running soft filtering..."

    bcftools view --threads 8 -m2 -Ou "$1/gatk.snp.indel.mnv.vcf.gz" | \
    bcftools filter --threads 8 -S . -e "FORMAT/DP<$minDP || FORMAT/DP>$maxDP" -Ou | \
    bcftools filter -s "(1)DRAGENHardQUAL" -e 'QUAL < 10.4139' -Ou | \
    bcftools filter -m + -s "(2)MISSING" -e 'F_MISSING > 0.2' -Ou | \
    bcftools filter -m + -s "(3)INVARIABLEREF" -e 'AC == 0' -Ou | \
    bcftools filter -m + -s "(4)MINMAF" -e "AF<$minmaf || AF>$minmafi" -Ou | \
    bcftools filter -m + -s "(5)MINMQ" -e "MQ < 40.0" -Ou | \
    bcftools filter -m + -s "(6)MINQD2" -e "QD < 2.0" -Ou | \
    bcftools filter -m + -s "(7)STRANDBIAS" -e "FS > 20" -Ou | \
    bcftools filter -m + -s "(8)MAXSOR3" -e "SOR > 3.0" -Ou | \
    bcftools filter -m + -s "(9)MQRankSum" -e "MQRankSum < -2.0 || MQRankSum > 4.0" -Ou | \
    bcftools filter -m + -s "(10)ReadPosRankSum" -e "ReadPosRankSum < -3.0 || ReadPosRankSum > 3.0" \
    -Oz -o "$SOFT_TMP"
fi

# Apply mask if provided
if [[ -n "${8:-}" && -f "$8" ]]; then
    echo "Applying mask..."
    bcftools filter --threads 8 -m + -s "(11)MASK" \
        --mask-file "$8" \
        -Oz -o "$SOFT_FINAL" "$SOFT_TMP"
else
    mv "$SOFT_TMP" "$SOFT_FINAL"
fi

# =====================================================
# 4. HARD FILTER
# =====================================================

if [[ ! -f "$HARD_FINAL" ]]; then
    echo "Creating hard filtered VCF..."

    bcftools filter --threads 8 \
        -i 'FILTER="PASS"' \
        -Oz -o "$HARD_FINAL" "$SOFT_FINAL"
fi

bcftools index -f "$SOFT_FINAL"
bcftools index -f "$HARD_FINAL"

# =====================================================
# 5. FILTER DISTRIBUTION
# =====================================================

echo "Calculating filter distribution..."

bcftools query -f '%CHROM\t%FILTER\n' "$SOFT_FINAL" | \
sort | uniq -c | \
awk '{print $2"\t"$3"\t"$1}' \
> "$1/gatk.snp.indel.mnv.filtered_soft.vcf.FILT_DIST"

# =====================================================
# 6. STATS
# =====================================================

echo "Running bcftools stats..."

bcftools stats \
    -v -d 5,100,5 \
    -S "$3" \
    "$HARD_FINAL" \
    > "$1/gatk.snp.indel.mnv.filtered.stats"

grep '^PSC' "$1/gatk.snp.indel.mnv.filtered.stats" | tail -n +2 \
> "$1/gatk.snp.indel.mnv.filtered.per-indiv.stats"

# =====================================================
# DONE
# =====================================================

end=$(date +%s)
echo "SUCCESS"
echo "Runtime: $((end-start)) sec"
