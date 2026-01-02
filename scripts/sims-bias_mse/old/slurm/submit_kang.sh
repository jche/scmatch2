#!/bin/bash
#SBATCH -J kang_sim
#SBATCH -o logs/kang_%A_%a.out
#SBATCH -e logs/kang_%A_%a.err
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=1-100

# Setup
project_dir="/homes2/xmeng/scmatch2/"
sim_type="kang"

# Load R module
module load R/4.3.2

# Create logs directory
mkdir -p logs

# Print job info
echo "========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Simulation: ${sim_type}"
echo "Starting at: $(date)"
echo "========================================="

# Change to project directory
cd ${project_dir}

# Run R script
Rscript scripts/sims-bias_mse/canonical_run_single_iteration.R ${sim_type} ${SLURM_ARRAY_TASK_ID}

echo "Finished at: $(date)"
