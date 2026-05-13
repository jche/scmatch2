#!/bin/bash
#SBATCH --job-name=var-multi
#SBATCH --array=1-998
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/var-multi_%A_%a.out
#SBATCH --error=logs/var-multi_%A_%a.err

mkdir -p logs

module load R/4.3.0   # adjust to your cluster's module name

Rscript scripts/sims-variance-multi/1_run_iter.R "${SLURM_ARRAY_TASK_ID}"
