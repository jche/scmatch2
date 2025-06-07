#!/bin/bash
#SBATCH -J parallel_sim        # Job name
#SBATCH -p sapphire            # Partition (queue) name (adjust if needed)
#SBATCH -o logs/parallel_sim_%A_%a.out  # Standard output log file (%A = Job ID, %a = Array Task ID)
#SBATCH -e logs/parallel_sim_%A_%a.err  # Standard error log file (%A = Job ID, %a = Array Task ID)
#SBATCH -c 1                   # Number of CPU cores per task
#SBATCH --mem=4G              # Memory per task (e.g., 4GB); adjust as needed
#SBATCH -t 01:00:00            # Max runtime per task (HH:MM:SS); adjust as needed
#SBATCH --array=1-500          # Job array: run tasks 1 through 500 (adjust range as needed)

# --- Environment Setup ---
my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18 # Adjust R library path if different
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif" # Adjust Singularity image path if different
project_dir="/n/holylabs/pillai_lab/Users/xmeng1/scmatch2/"
r_script_path="${project_dir}/scripts/inference-scripts/1.5_parallel_sim_inference.R" # Path to the new R script

# Create logs directory if it doesn't exist
mkdir -p "${project_dir}/scripts/inference-scripts/logs"

# Check if R script exists
if [ ! -f "$r_script_path" ]; then
  echo "Error: R script not found at $r_script_path"
  exit 1
fi


singularity_command="singularity exec \
  --cleanenv \
  --env R_LIBS_USER=${my_packages} \
  --bind ${project_dir}:${project_dir} \
  ${rstudio_singularity_image}"


echo "Starting R script for iteration ${SLURM_ARRAY_TASK_ID}"
cd ${project_dir}

$singularity_command Rscript ${r_script_path} ${SLURM_ARRAY_TASK_ID}

echo "Finished R script for iteration ${SLURM_ARRAY_TASK_ID}"

