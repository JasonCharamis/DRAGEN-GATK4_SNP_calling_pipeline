#!/bin/bash
#SBATCH --job-name=GATK4genotype

scaff=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1/scaffolds.txt)

start=$(date +%s)
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo ${scaff}

dc=$(cat $1/conda.txt | awk '{print $1}')
ec=$(cat $1/conda.txt | awk '{print $2}')

source ${dc}/conda.sh
conda activate ${ec}

# GenotypeGVCFs
if [[ -f "$2/scaff_vcfs/${scaff}.gatk.vcf.gz" ]]; then
	echo "$2/scaff_vcfs/${scaff}.gatk.vcf.gz exists. Skipping joint genotyping."
else
	echo "Joint genotyping all samples..."
	gatk --java-options "-Xmx7g" GenotypeGVCFs \
		--reference $3 \
		--variant gendb://$2/genomicDB/${scaff}.gvcf.db \
		--output $2/scaff_vcfs/${scaff}.gatk.vcf.gz \
		--tmp-dir $1
fi

if [[ -f $1/${scaff}.vcf.rehead ]]; then
	echo "$1/${scaff}.vcf.rehead exists. Skipping rehead."
else
	zcat $2/scaff_vcfs/${scaff}.gatk.vcf.gz | grep -m 1 "#CHR" | cut -f 10- | tr "\t" "\n" | sed -E "s/\_$scaff//" >$1/${scaff}.vcf.rehead
fi

(exit) && echo success
end=$(date +%s)
runtime=$((end - start))
echo "RUNTIME: $runtime"
