#!/bin/bash
set -e

#!/bin/bash
#SBATCH -A account_name   # account name
#SBATCH -D working_directory  # working directory
#SBATCH --job-name=variant_calling_pipeline_DRAGEN_GATK4_BP
#SBATCH --time=20:00:00
#SBATCH --output=00_DRAGEN_GATK_pipeline_%j.out
#SBATCH --error=00_DRAGEN_GATK_pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --exclude=comp[18-24]

###############
# DESCRIPTION #
###############
# This pipeline performs whole-genome variant calling (SNPs, MNVs, indels) of whole genome resequencing data following the DRAGEN-GATK best practices protocol (https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939)
#  with the exception of extra preprocessing utilizing HTStream. Preprocessing largely follows the recommendations of the UC Davis Bioinformatics Core
#  (https://ucdavis-bioinformatics-training.github.io/2020-mRNA_Seq_Workshop/data_reduction/01-preproc_htstream_mm) with the addition of filtering reads with short (< 100bp) insert lengths, a necessary#  step for avoidance of segmenation fault errors for DRAGMAP v1.2.1. This pipeline works with samples sequenced on individual lanes or samples sequenced across multiple lanes. Restriction of SNP-calling to specific genomic intervals as well as masking of unwanted regions is also possible by specifying whitelist or blacklist files (see below). This pipeline was developed using the job scheduling
#  manager SLURM on an HPC cluster.

################
# DEPENDENCIES #
################
# HTStream v1.4.1       (https://s4hts.github.io/HTStream/#hts_QWindowTrim; conda install -c bioconda htstream)
# DRAGMAP v1.3.0 (https://github.com/Illumina/DRAGMAP; conda install -c bioconda dragmap=1.2.1)
# picard v2.20.4        (https://github.com/broadinstitute/picard; conda install -c bioconda picard)
# GATK4 v4.6.2.0          (https://github.com/broadinstitute/gatk; conda install -c bioconda gatk4)
# samtools v1.23       (https://github.com/samtools/samtools; conda install -c bioconda samtools)
# bedtools v2.31.1       (https://github.com/arq5x/bedtools2; conda install -c bioconda bedtools)
# bcftools v1.23      (https://github.com/samtools/bcftools; conda install -c bioconda bcftools)
## NOTE: All packages are loaded through conda

##################
# REQUIRED FILES #
##################
# raw fastq files
# reference genome fasta
# sample list(s): file containing sample prefixes (one per line); if reads need to be merged across lanes create two separate prefix lists - one with and one without lane designation
# read group information file: tab delimited file containing the following read group information for each (unmerged) sample (one per line) in the following order - ID SM LB PL PU
## ID = sample read group ID (e.g. sample prefix including lane designation; identical to SM if sequenced on a single lane)
## SM = sample name (excluding lane designation)
## LB = DNA prep library identifyier (only imporatant if multiple libraries were sequenced)
## PU = platform unit - code(s) for the sequencing unit/machine and lane used; if the same sequencing unit was used for all libraries just put lane designation (e.g. L001)
# OPTIONAL: bed specifying genomic interval(s) to include OR bed file specifying genomic interval(s) to exclude/mask

#################
# SET VARIABLES #
#################
dir="" # path to root directory for all input and output files

dirconda=""         # path to directory containing conda.sh
envconda="genomics" # conda environment containing DRAGMAP & picard (& fastq-filter?)

ref="${dir}/reference/GCF_000001635.27_GRCm39_genomic.fna" # path to reference genome
mtDNA="NC_005089.1"                                        # name of mtDNA scaffold if known
include=""                                                 # OPTIONAL bed file with intervals to detect SNPs within; for whole genome leave blank
exclude=""                                                 # OPTIONAL path to masking bed file for SNP calling

raw="${dir}/fastq"                             # directory containing raw reads
samples="${dir}/fastq/samples_list.txt"        # path to file containing fastq sample prefixes (including individual lane designations if applicable)
samples_merged="${dir}/fastq/samples_list.txt" # path to file containing merged sample prefixes (i.e. excluding lane designations); if samples sequenced on single lane set to same file path as ${samples}
rg="${dir}/fastq/rg_info.tsv"                  # read group information file; see README.md for details
bin_dir="${dir}/bin"

## SET THE FOLLOWING TO 'T' OR 'F' to delete temporary and intermediate files##
rmtemp="T" # set to T to delete the directory containing temporary files (i.e. ${dir}/temp)
rmhts="F"  # set to T to delete cleaned FASTQs (i.e. ${dir}/01_HTS_PreProc/*fastq)
rmbam="F"  # set to T to delete bams (i.e. ${dir}/02_AlignDRAGMAP/*.bam)
rmivcf="F" # set to T to delete gvcfs and individual scaffold vcfs (i.e. ${dir}/04_GATKvcfs/indiv_scaff_gvcfs AND ${dir}/04_GATKvcfs/scaff_vcfs)

## IF READS NEED TO BE MERGED ACROSS LANES ##
merge="" # set to any string if fastqs need to be merged across multiple lanes, if not leave empty
lane=""  # shared lane designation prefix (e.g. "_L00" for sample_prefix_L001, sample_prefix_L002, etc.)

## NOTE: time limits and/or memory allocation may need to be adjusted depending on the size of your data

# Optimized for a node with 40 cores / RAM

# Load modules / conda
source /miniconda3/etc/profile.d/conda.sh
conda init
conda activate ${envconda}

# -------------------------
# Pipeline directories
# -------------------------
temp="${dir}/temp"
std="${dir}/slurmout"
preproc="${dir}/01_HTS_PreProc"
align="${dir}/02_AlignDRAGMAP"
stats="${dir}/03_AlignStats"
vcfs="${dir}/04_GATKvcfs"

mkdir -p ${std}/{01_htstream,02_DRAGMAP,03_merge,04_stats,05_hapcall,06_genotype}
mkdir -p ${dir}/temp ${preproc} ${align} ${stats} ${vcfs}/dragstr_model ${vcfs}/indiv_scaff_gvcfs ${vcfs}/scaff_vcfs ${vcfs}/genomicDB

# Create file with directory and conda environment
echo -e "${dirconda}\t${envconda}" >"${temp}/conda.txt"

# Index genome using samtools
samtools faidx ${ref}

# Create scaffolds list
if [ ${include} ]; then
	cat ${include} | cut -f1 | sort | uniq >${dir}/reference/scaffolds.list
else
	cat ${ref}.fai | cut -f1 >${dir}/reference/scaffolds.list
fi

# Create file with all individual × scaffold combinations
cat $samples_merged | while read indiv; do
	cat ${dir}/reference/scaffolds.list | while read scaff; do
		echo -e "${indiv}\t${scaff}"
	done
done >${temp}/indiv_scaff.tsv

# Add library info from read group file
cat ${rg} | cut -f2-3 | sort | uniq | while read line; do
	i=$(echo $line | awk '{print $1}')
	l=$(echo $line | awk '{print $2}')
	cat ${temp}/indiv_scaff.tsv | grep -w ${i} | awk -v lib=$l '{print $0"\t"lib}'
done >${temp}/temp.tsv && mv ${temp}/temp.tsv ${temp}/indiv_scaff.tsv

# Count samples / scaffolds
total=$(cat $samples | wc -l)
nindiv=$(cat $samples_merged | wc -l)
nscaffs=$(cat ${dir}/reference/scaffolds.list | cut -f2 | sort | uniq | wc -l)
ind_scaf=$((${nindiv} * ${nscaffs}))

# -------------------------
# CLEAN & ALIGN
# -------------------------
# 1. Preprocessing
htsid=$(sbatch --parsable -t 5-00:00:00 \
	--array=1-${total}%50 \
	--nodes=1 \
	--ntasks=1 \
	--cpus-per-task=7 \
	--mem=23G \
	--output=${std}/01_htstream/01_htstream_%A_%a.out \
	--error=${std}/01_htstream/01_htstream_%A_%a.err \
	${bin_dir}/HTS_preproc.slurm ${samples} ${preproc} ${raw} ${temp})

# 2. Hash Table
htid=$(sbatch --parsable -t 10-00:00:00 \
	--nodes=1 \
	--ntasks=1 \
	--cpus-per-task=36 \
	--mem=112G \
	--output=${std}/dragmap_hash_%A.out \
	--error=${std}/dragmap_hash_%A.err \
	${bin_dir}/hashDRAGMAP.slurm ${temp} ${ref})

# 3. DRAGMAP Alignment
dmid=$(sbatch --parsable -t 5-00:00:00 \
	--array=1-${total}%9 \
	--ntasks=1 \
	--cpus-per-task=13 \
	--mem=38G \
	--dependency=afterok:${htsid}:${htid} \
	--output=${std}/02_DRAGMAP/02_DRAGMAP_%A_%a.out \
	--error=${std}/02_DRAGMAP/02_DRAGMAP_%A_%a.err \
	${bin_dir}/alignDRAGMAP.slurm ${samples} ${align} ${ref} ${temp} ${preproc} ${rg} ${lane} ${rmhts})

# 4. Genome windows
gwid=$(sbatch -t 01-00:00:00 \
	--mem=1G \
	--output=${std}/genome_win_%A.out \
	--error=${std}/genome_win_%A.err \
	${bin_dir}/genome_wins.slurm ${ref} ${temp} | sed 's/Submitted batch job //')

# 5. Merge BAMs (if samples were sequenced in different lanes)
if [ ${merge} ]; then
	mgid=$(sbatch -t 3-00:00:00 \
		--array=1-${nindiv}%5 \
		--output=${std}/03_merge/03_merge_%A_%a.out \
		--error=${std}/03_merge/03_merge_%A_%a.err \
		--dependency=afterok:${dmid} \
		${bin_dir}/samtools_merge.slurm ${samples_merged} ${align} ${lane} ${temp} | sed 's/Submitted batch job //')

	dep_asid="${mgid}:${gwid}"
else
	dep_asid="${dmid}:${gwid}"
fi

# 6. Alignment statistics
asid=$(sbatch -t 3-00:00:00 \
	--array=1-${nindiv}%10 \
	--output=${std}/04_stats/04_stats_%A_%a.out \
	--error=${std}/04_stats/04_stats_%A_%a.err \
	--dependency=afterok:${dep_asid} \
	${bin_dir}/align_stats.slurm ${samples_merged} ${stats} ${ref} ${align} ${temp} | sed 's/Submitted batch job //')

# -------------------------
# Variant Calling - in its current form SNPs, MNVs (multi-nucleotide variants) and indels are selected, but you can change that in bin/vcf_scaff_to_snp.vcf.slurm
# -------------------------
# 7. STR table
stid=$(sbatch -t 2-00:00:00 \
	--output=${std}/str_table_%A.out \
	--error=${std}/str_table_%A.err \
	${bin_dir}/STRtable.slurm ${ref} ${temp} ${include} | sed 's/Submitted batch job //')

# 8. Haplotype calling (5 runs per node with 8 cpus each, in 10 nodes)
htcid=$(sbatch -t 10-00:00:00 \
	--array=1-${ind_scaf}%50 \
	--nodes=1 \
	--ntasks=1 \
	--cpus-per-task=8 \
	--mem=23G \
	--error=${std}/05_hapcall/05_hapcall_%A_%a.err \
	--output=${std}/05_hapcall/05_hapcall_%A_%a.out \
	--dependency=afterok:${stid}:${asid} \
	${bin_dir}/bam_to_gvcf.slurm ${temp} ${ref} ${align} ${vcfs} ${rmbam} ${include} | sed 's/Submitted batch job //')

# Splitted GenomicsDB import and GenotypeVCFs because of different computing needs; GenomicsDB import will use many threads for 2 jobs per node at a time (for GATK4's reader threads)
# while GenotypeVCFs which has a lot of I/O processes (and therefore multiple threads wouldn't translate into smaller running time)

# 9. GenomicsDB import
gdbid=$(sbatch -t 10-00:00:00 \
	--array=1-${nscaffs}%16 \
	--nodes=1 \
	--ntasks=1 \
	--mem=55G \
	--cpus-per-task=20 \
	--error=${std}/06_genotype/06_genomicsdb_import_%A_%a.err \
	--output=${std}/06_genotype/06_genomicsdb_import_%A_%a.out \
	--dependency=afterok:${htcid} \
	${bin_dir}/db_import.slurm ${temp} ${vcfs} ${ref} | sed 's/Submitted batch job //')

# 10. Joint genotyping
jgtid=$(sbatch -t 10-00:00:00 \
	--array=1-${nscaffs}%50 \
	--nodes=1 \
	--ntasks=1 \
	--mem=5G \
	--cpus-per-task=2 \
	--error=${std}/06_genotype/06_joint_genotyping_%A_%a.err \
	--output=${std}/06_genotype/06_joint_genotyping_%A_%a.out \
	--dependency=afterok:${gdbid} \
	${bin_dir}/genotype_vcfs.slurm ${temp} ${vcfs} ${ref} | sed 's/Submitted batch job //')

# 11. Concatenate & filter variants
csid=$(sbatch -t 5-00:00:00 \
	--error=${std}/07_concat_vcfs_filter_SNPs_%A.err \
	--output=${std}/07_concat_vcfs_filter_SNPs_%A.out \
	--dependency=afterok:${jgtid} \
	${bin_dir}/vcf_scaff_to_snp.vcf.slurm ${vcfs} ${temp} ${samples_merged} ${mtDNA} ${nindiv} ${rmivcf} ${exclude} | sed 's/Submitted batch job //')

# 12. Job report
sucid=$(sbatch -t 01-00:00:00 \
	--error=${std}/JOB_REPORT_%A.err \
	--output=${std}/JOB_REPORT_%A.out \
	--dependency=after:${csid} \
	${bin_dir}/job_report.slurm ${std} ${temp} ${rmtemp} | sed 's/Submitted batch job //')

# Define your permanent home results directory
home_results="/home/gsentis/jason/gender_bias/variant_calling_results_$(date +%F)"

echo -e "${htsid}\n${htid}\n${dmid}\n${gwid}\n${asid}\n${stid}\n${htcid}\n${jgtid}\n${csid}\n${mgid}" >${temp}/job_IDs.txt

(exit) && echo success
