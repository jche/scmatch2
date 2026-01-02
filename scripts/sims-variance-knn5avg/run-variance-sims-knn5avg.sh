#!/bin/bash
#SBATCH -J csm_infer_knn5avg
#SBATCH -o logs/csm_infer_knn5avg_%A_%a.out
#SBATCH -e logs/csm_infer_knn5avg_%A_%a.err
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --array=1-1000

project_dir="/homes2/xmeng/scmatch2"
r_script_path="${project_dir}/scripts/sims-variance-knn5avg/1_parallel_sim_inference.R"

module load R/4.3.2

mkdir -p "${project_dir}/logs"

if [ ! -f "$r_script_path" ]; then
  echo "Error: R script not found at $r_script_path"
  exit 1
fi

echo "Starting R script for iteration ${SLURM_ARRAY_TASK_ID}"
cd "${project_dir}"

Rscript "${r_script_path}" "${SLURM_ARRAY_TASK_ID}"

echo "Finished iteration ${SLURM_ARRAY_TASK_ID}"
