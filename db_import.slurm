#!/bin/bash
#SBATCH --job-name=GenomicsDBImport

scaff=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1/scaffolds.txt)

start=$(date +%s)
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo ${scaff}

dc=$(cat $1/conda.txt | awk '{print $1}')
ec=$(cat $1/conda.txt | awk '{print $2}')

source ${dc}/conda.sh
conda activate ${ec}

ls -d $2/indiv_scaff_gvcfs/*${scaff}.rb.g.vcf.gz | sed -E "s/(.+\/)(.+)\_${scaff}\.rb.g.vcf.gz/\2\t\1\2\_${scaff}.rb.g.vcf.gz/" >$1/samples_${scaff}.gvcf.txt

# GenomicsDBImport
echo "Importing single-sample gvcfs into database..."
gatk --java-options "-Xmx45g -Xms45g" GenomicsDBImport \
	--genomicsdb-workspace-path $2/genomicDB/${scaff}.gvcf.db \
	--intervals $1/${scaff}.bed \
	--sample-name-map $1/samples_${scaff}.gvcf.txt \
	--reader-threads 10 \
	--batch-size 50 \
	--use-jdk-deflater \
	--use-jdk-inflater \
	--tmp-dir $1

end=$(date +%s)
runtime=$((end - start))
echo "RUNTIME: $runtime"
