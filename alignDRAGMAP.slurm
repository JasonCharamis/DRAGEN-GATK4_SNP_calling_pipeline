#!/bin/bash

#SBATCH --job-name=alignDRAGMAP

# --- 1. Variables and Setup ---
sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" "$1")
#ind=$(echo "${sample}" | sed -E "s/$7.+//")
ref=$3
refdir=$(dirname "${ref}")/dragen_hash

# Extract metadata with safety fallbacks
metadata_row=$(grep -w "${sample}" "$6")
sm=$(echo "${metadata_row}" | awk '{print $2}')
lib=$(echo "${metadata_row}" | awk '{print $3}')
pu=$(echo "${metadata_row}" | awk '{print $4}')

# Environment setup
dc=$(cat "$4/conda.txt" | awk '{print $1}')
ec=$(cat "$4/conda.txt" | awk '{print $2}')

echo "Host: $HOSTNAME"
echo "Sample: ${sample}"
echo "Array ID: $SLURM_ARRAY_TASK_ID"

[[ -d "$2/${sample}" ]] || mkdir -p "$2/${sample}"

source "${dc}/conda.sh"
conda activate "${ec}"
export TMPDIR=$4

# If the RG BAM exists, the whole job is done. If the sorted BAM exists, skip alignment.
if [[ -f "$2/${sample}/${sample}.sorted.rg.md.bam" || \
      -f "$2/${sample}/${sample}.sorted.rg.bam"    || \
      -f "$2/${sample}/${sample}.sorted.bam" ]]; then
    echo "$2/${sample}/${sample}.sorted.bam already exists. Skipping alignment."
elif [[ ! -f "$2/${sample}/${sample}.sorted.bam" ]]; then
    echo "Aligning reads..."
    dragen-os -v \
        -r "${refdir}" \
        -1 "$5/${sample}/${sample}_preproc_R1.fastq.gz" \
        -2 "$5/${sample}/${sample}_preproc_R2.fastq.gz" \
        --num-threads $SLURM_CPUS_PER_TASK | \
        samtools view -bh - | \
        samtools sort -@ 2 -O bam -o "$2/${sample}/${sample}.sorted.bam"
else
    echo "Found existing sorted.bam. Skipping alignment, proceeding to Read Groups."
fi

if [[ -f "$2/${sample}/${sample}.sorted.rg.md.bam" || \
      -f "$2/${sample}/${sample}.sorted.rg.bam" ]]; then

    echo "$2/${sample}/${sample}.sorted.rg.bam already exists. Skipping adding read groups."

else
    echo "Replacing read groups..."
    
    # If metadata variables are empty, assign defaults so GATK doesn't crash
    [[ -z "${sm}" ]]  && sm="${sample}"
    [[ -z "${lib}" ]] && lib="lib1"
    [[ -z "${pu}" ]]  && pu="run1"

    gatk AddOrReplaceReadGroups \
        -I "$2/${sample}/${sample}.sorted.bam" \
        -O "$2/${sample}/${sample}.sorted.rg.bam" \
        -PL Illumina \
        -LB "${lib}" \
        -SM "${sm}" \
        -PU "${pu}" && \
	rm "$2/${sample}/${sample}.sorted.bam"
fi

if [[ -f "$2/${sample}/${sample}.sorted.rg.md.bam" ]]; then
    echo "MarkDuplicates has been performed. Skipping this step."
else
    gatk MarkDuplicates \
	 -I $2/${sample}/${sample}.sorted.rg.bam \
	 -O $2/${sample}/${sample}.sorted.rg.md.bam \
	 -M $2/${sample}/${sample}.sorted.rg.md.metric.txt && \
	rm $2/${sample}/${sample}.sorted.rg.bam
fi

echo "Job for ${sample} complete."
