#!/bin/bash

#SBATCH --job-name=align_stats

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
[[ -d $2/${sample} ]] || mkdir -p $2/${sample}
ref=$3
refpre=$(echo ${ref%.*})

start=$(date +%s)
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${sample}"

# Load samtools and bedtools from conda
dc=$(cat $5/conda.txt | awk '{print $1}')
ec=$(cat $5/conda.txt | awk '{print $2}')

source ${dc}/conda.sh
conda activate ${ec}

if [[ -f $4/${sample}/${sample}.sorted.rg.md.bam.bai ]]; then
	echo "$4/${sample}/${sample}.sorted.rg.md.bam exists. Skipping indexing it."
else
	echo "indexing bam..."
	samtools index $4/${sample}/${sample}.sorted.rg.md.bam
fi

echo "Calculating samtools stats:"

if [[ -f "$2/${sample}/${sample}.dragmap_flagstat.txt" ]]; then
	echo "$2/${sample}/${sample}.dragmap_flagstat.txt exists. Skipping flagstat."
else
	echo "  flagstat..."
	samtools flagstat $4/${sample}/${sample}.sorted.rg.md.bam >$2/${sample}/${sample}.dragmap_flagstat.txt
fi

if [[ -f "$2/${sample}/${sample}.samtools_stats.txt" ]]; then
	echo "$2/${sample}/${sample}.samtools_stats.txt exists. DRAGMAP stats skipped."
else
	echo "  stats..."
	samtools stats -c 1,10000,100 $4/${sample}/${sample}.sorted.rg.md.bam >$2/${sample}/${sample}.samtools_stats.txt
fi

if [[ -f "$2/${sample}/${sample}.dragmap_IDstats.tsv" ]]; then
	echo "$2/${sample}/${sample}.dragmap_IDstats.tsv exists. Skipping idxstats."
else
	echo "  idxstats..."
	samtools idxstats $4/${sample}/${sample}.sorted.rg.md.bam >$2/${sample}/${sample}.dragmap_IDstats.tsv
fi

if [[ -f "$2/${sample}/${sample}.dragmap_depth.txt" ]]; then
	echo "$2/${sample}/${sample}.dragmap_depth.txt exists. Skipping DRAGMAP depth estimation."
else
	echo "  depth..."
	samtools depth $4/${sample}/${sample}.sorted.rg.md.bam | awk '{sum+=$3; cnt++}END{print "no_positions\tmean_depth"; print sum"\t"sum/cnt}' >$2/${sample}/${sample}.dragmap_depth.txt
fi

if [[ -f "$2/${sample}/${sample}.dragmap_coverage.tsv" ]]; then
	echo "$2/${sample}/${sample}.dragmap_coverage.tsv exists. Skipping DRAGMAP coverage estimation."
else
	echo "  coverage..."
	samtools coverage $4/${sample}/${sample}.sorted.rg.md.bam >$2/${sample}/${sample}.dragmap_coverage.tsv
fi

if [[ -f "$2/${sample}/${sample}.dragmap.gatk.stat.txt" ]]; then
	echo "$2/${sample}/${sample}.dragmap.gatk.stat.txt exists. Skipping GATK metrics estimation."
else
	echo "Calculating GATK metrics..."
	echo "checking mapping quality"
	gatk CollectAlignmentSummaryMetrics \
		-R $3 \
		-I $4/${sample}/${sample}.sorted.rg.md.bam \
		-O $2/${sample}/${sample}.dragmap.gatk.stat.txt
fi

if [[ -f "$2/${sample}/${sample}.dragmap_50kbwin_meancov.tsv" ]]; then
	echo "$2/${sample}/${sample}.dragmap_50kbwin_meancov.tsv exists. Skipping DRAGMAP 50kb window mean coverage estimation."
else
	echo "Calculating coverage for 50Kb windows..."
	echo -e "scaff\twin_start\twin_end\tmean_depth" >$2/${sample}/${sample}.dragmap_50kbwin_meancov.tsv
	bedtools coverage -mean -sorted -a ${refpre}_50kb_win.bed -b $4/${sample}/${sample}.sorted.rg.md.bam >>$2/${sample}/${sample}.dragmap_50kbwin_meancov.tsv
fi

if [[ -f "$2/${sample}/${sample}.coverage_histogram.txt" ]]; then
	echo "$2/${sample}/${sample}.coverage_histogram.txt exists. Skipping coverage estimation using bedtools genomecov."
else
	# 1. Filter out duplicates (-F 1024)
	# 2. Compute coverage with bedtools
	# 3. Create a BedGraph (-bg) for IGV and a histogram (-max) for statistics
	# We use 'tee' to split the stream.
	# One side goes to a 'process substitution' for the BedGraph
	# The other side continues to the Histogram.
	samtools view -u -F 1024 ${sample}/${sample}.sorted.rg.md.bam |
		tee >(bedtools genomecov -ibam stdin -bg >${sample}.no_dups.bedgraph) |
		bedtools genomecov -ibam stdin -max 100 >${sample}.coverage_histogram.txt
fi

(exit) && echo success
end=$(date +%s)
runtime=$((end - start))
echo "RUNTIME: $runtime seconds"
