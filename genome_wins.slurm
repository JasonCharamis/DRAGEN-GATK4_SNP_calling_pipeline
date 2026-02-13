#!/bin/bash

#SBATCH --job-name=genome_win

start=`date +%s`
echo $HOSTNAME

ref=$1
refpre=$(echo ${ref%.*})

dc=$(cat $2/conda.txt | awk '{print $1}')
ec=$(cat $2/conda.txt | awk '{print $2}')

source ${dc}/conda.sh
#conda init
conda activate ${ec}

if [[ -f "${refpre}_50kb_win.bed" ]]; then
    echo "${refpre}_50kb_win.bed exists. Skipping this step."
else 
    echo "breaking genome into 50kb windows..."
    cat $1.fai | awk '{print $1"\t"0"\t"$2-1}' > $refpre.bed
    bedtools makewindows -b $refpre.bed -w 50000 -s 45000 > ${refpre}_50kb_win.bed
fi
    
(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime seconds"
