#!/bin/bash
#SBATCH -J acic_sim           # Job name
#SBATCH -o logs/acic_%A_%a.out  # Standard output log
#SBATCH -e logs/acic_%A_%a.err  # Standard error log
#SBATCH -c 4                   # Number of CPU cores per task
#SBATCH --mem=8G               # Memory per task
#SBATCH -t 02:00:00            # Max runtime (2 hours per iteration)
#SBATCH --array=1-100          # Run 100 iterations

# Setup
project_dir="/homes2/xmeng/scmatch2/"
sim_type="acic"

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
