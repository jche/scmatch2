#!/bin/bash
#SBATCH -J test_run        # Job name
#SBATCH -o logs/test_run_%A_%a.out  # Standard output log file (%A = Job ID, %a = Array Task ID)
#SBATCH -e logs/test_run_%A_%a.err  # Standard error log file (%A = Job ID, %a = Array Task ID)
#SBATCH -c 60                   # Number of CPU cores per task
#SBATCH --mem=24G              # Memory per task (e.g., 4GB); adjust as needed
#SBATCH -t 01:00:00            # Max runtime per task (HH:MM:SS); adjust as needed
#SBATCH --array=1-500          # Job array: run tasks 1 through 500 (adjust range as needed)

project_dir="/homes2/xmeng/scmatch2/"
r_script_path="${project_dir}/scripts/inference-scripts/1_parallel_sim_inference.R"

module load R/4.3.2

# Set a reference to /homes2/xmeng/R/x86_64-pc-linux-gnu-library/4.4
# There is another package installation path on kraken -- how to find that
# export R_LIBS_USER=/homes2/xmeng/R/x86_64-pc-linux-gnu-library/4.4

# Create logs directory if it doesn't exist
mkdir -p "${project_dir}/scripts/inference-scripts/logs"


# Check if R script exists
if [ ! -f "$r_script_path" ]; then
  echo "Error: R script not found at $r_script_path"
  exit 1
fi


echo "Starting R script for iteration ${SLURM_ARRAY_TASK_ID}"
cd ${project_dir}

Rscript ${r_script_path} ${SLURM_ARRAY_TASK_ID}

echo "Finished R script for iteration ${SLURM_ARRAY_TASK_ID}"

