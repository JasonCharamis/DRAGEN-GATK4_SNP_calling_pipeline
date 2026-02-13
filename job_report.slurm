#!/bin/bash
#SBATCH --job-name=JOB_REPORT
#SBATCH --output=job_report.out
#SBATCH --mem=1G

echo "Pipeline finished at $(date)"
ls -R $1 $2
