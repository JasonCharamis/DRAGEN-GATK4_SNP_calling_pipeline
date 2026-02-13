#!/bin/bash

#SBATCH --job-name=merge

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`

start=`date +%s`
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${sample}"

dc=$(cat $4/conda.txt | awk '{print $1}')
ec=$(cat $4/conda.txt | awk '{print $2}')

source ${dc}/conda.sh
conda activate ${ec}

if [[ -f "$2/${sample}/${sample}.sorted.rg.bam" ]]; then
    echo "$2/${sample}/${sample}.sorted.rg.bam exists. Skipping merging sample reads across lanes."
else
    echo "Merging sample reads across lanes..."
    samtools merge -f -o $2/${sample}/${sample}.sorted.rg.bam $2/${sample}/${sample}*$3*.sorted.rg.bam
    rm $2/${sample}/*$3*.bam
fi

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime seconds"

