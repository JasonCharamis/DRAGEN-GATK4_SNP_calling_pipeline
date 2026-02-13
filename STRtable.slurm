#!/bin/bash

#SBATCH --job-name=str_table

ref=$1
refdir=$(echo ${ref%/*})
refpre=$(echo ${ref%.*})
dc=$(cat $2/conda.txt | awk '{print $1}')
ec=$(cat $2/conda.txt | awk '{print $2}')

start=`date +%s`
echo $HOSTNAME

source ${dc}/conda.sh
conda activate ${ec}

if [[ -f "$refpre.dict" ]]; then
    echo "$refpre.dict exists. Skipping creation of sequence dictionary".
else
    echo "Indexing genome (gatk)..."
    gatk CreateSequenceDictionary -R $1 -O $refpre.dict
fi

if [[ -f "$refdir/str_table.tsv" ]]; then
    echo "$refdir/str_table.tsv exists. Skipping STR table composing".
else
    echo "Composing STR table..."
    if [ $3 ]
    then
	gatk ComposeSTRTableFile \
	     -R $1 \
	     -O $refdir/str_table.tsv \
	     -L $3
    else
	gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' ComposeSTRTableFile \
	     -R $1 \
	     -O $refdir/str_table.tsv
    fi
fi

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
