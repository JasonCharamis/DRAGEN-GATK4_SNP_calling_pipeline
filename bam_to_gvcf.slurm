#!/bin/bash

#SBATCH --job-name=GATK4hapcall

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" "$1/indiv_scaff.tsv")
ref=$2
refdir=${ref%/*}

samp=$(echo "$sample" | awk '{print $1}')
scaff=$(echo "$sample" | awk '{print $2}')

mkdir -p "$1/${samp}_${scaff}"

dc=$(awk '{print $1}' "$1/conda.txt")
ec=$(awk '{print $2}' "$1/conda.txt")

source "${dc}/conda.sh"
conda activate "${ec}"

start=$(date +%s)

echo "$HOSTNAME"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "$sample"

FINAL_GVCF="$4/indiv_scaff_gvcfs/${samp}_${scaff}.rb.g.vcf.gz"
RAW_GVCF="$4/indiv_scaff_gvcfs/${samp}_${scaff}.g.vcf.gz"
DRAGSTR_MODEL="$4/dragstr_model/${samp}_${scaff}.dragstr_model.txt"

# If final reblocked gVCF exists → job is done
if [[ -f "$FINAL_GVCF" ]]; then
	echo "Final reblocked gVCF exists. Job already complete."
	exit 0
fi

# Ensure DRAGSTR model exists
if [[ ! -f "$DRAGSTR_MODEL" ]]; then
	echo "Calibrating DRAGEN STR model..."

	gatk --java-options "-Xms16g -Xmx16g" CalibrateDragstrModel \
		-R "$2" \
		-I "$3/${samp}/${samp}.sorted.rg.md.bam" \
		-L "$refdir/scaffolds.list" \
		-str "$refdir/str_table.tsv" \
		-O "$DRAGSTR_MODEL"
else
	echo "DRAGSTR model exists. Skipping calibration."
fi

# Ensure raw gVCF exists
if [[ ! -f "$RAW_GVCF" ]]; then

	if [[ -n "$6" ]]; then
		echo "Calling haplotypes over specified intervals..."
		INTERVAL_OPTION="-L $6"
	else
		echo "Calling haplotypes over default intervals..."
		INTERVAL_OPTION="-L $1/scaffolds.list"
	fi

	gatk --java-options "-Xmx20g" HaplotypeCaller \
		-R "$2" \
		-I "$3/${samp}/*rg.md.bam" \
		-O "$RAW_GVCF" \
		$INTERVAL_OPTION \
		-ERC GVCF \
		--dragen-mode true \
		--native-pair-hmm-threads "$SLURM_CPUS_PER_TASK" \
		--dragstr-params-path "$DRAGSTR_MODEL" \
		--use-jdk-deflater \
		--use-jdk-inflater \
		--tmp-dir "$1/${samp}_${scaff}"

else
	echo "Raw gVCF exists. Skipping HaplotypeCaller."
fi

# Reblock gVCF
if [[ ! -f "$FINAL_GVCF" ]]; then
	echo "Reblocking gVCF..."

	gatk --java-options "-Xms16g -Xmx16g" ReblockGVCF \
		-R "$2" \
		-V "$RAW_GVCF" \
		-O "$FINAL_GVCF" \
		--use-jdk-deflater \
		--use-jdk-inflater
fi

# If final output exists remove BAM files if option selected
if [[ -f "$FINAL_GVCF" ]]; then
	echo "Final output confirmed. Cleaning up..."

	if [[ "$5" == "T" ]]; then
		rm -f "$3/${samp}/*rg.md.bam"
	fi

	rm -f "$RAW_GVCF"*
fi

(exit) && echo success
end=$(date +%s)
runtime=$((end - start))
echo "RUNTIME: $runtime seconds"
