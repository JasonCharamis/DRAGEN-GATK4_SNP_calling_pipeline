#!/bin/bash

#SBATCH --job-name=dragmap_hash

start=`date +%s`
echo $HOSTNAME

dc=$(cat $1/conda.txt | awk '{print $1}')
ec=$(cat $1/conda.txt | awk '{print $2}')

source ${dc}/conda.sh
conda activate ${ec}

ref=$2
refdir=$(dirname "${ref}")/dragen_hash
mkdir -p "${refdir}"

if [[ -f "${refdir}/hash_table.complete" ]]; then
    echo "DRAGMAP hash tables already built. Skipping."
else
    echo "Creating hash tables for DRAGMAP..."

    dragen-os \
        --build-hash-table true \
        --ht-reference "${ref}" \
        --output-directory "${refdir}" \
        --ht-num-threads $SLURM_CPUS_PER_TASK

    # Mark completion only if command succeeds
    if [[ $? -eq 0 ]]; then
        touch "${refdir}/hash_table.complete"
    fi
fi

end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"

(exit) && echo success
