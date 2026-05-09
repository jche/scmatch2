#!/bin/bash
#SBATCH -J csm_adaptive_k2_scm
#SBATCH -o logs/csm_adaptive_k2_scm_%A_%a.out
#SBATCH -e logs/csm_adaptive_k2_scm_%A_%a.err
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --array=1-1000

project_dir="/homes2/xmeng/scmatch2"
r_script_path="${project_dir}/scripts/sims-variance/1_parallel_sim_inference.R"

module load R/4.3.2

mkdir -p "${project_dir}/logs"

if [ ! -f "$r_script_path" ]; then
  echo "Error: R script not found at $r_script_path"
  exit 1
fi

echo "Starting iter=${SLURM_ARRAY_TASK_ID} [adaptive, k=2, scm]"
cd "${project_dir}"

Rscript "${r_script_path}" "${SLURM_ARRAY_TASK_ID}" "sims-variance-fully-adaptive(2nn)" "adaptive" 2 "scm"

echo "Finished iter=${SLURM_ARRAY_TASK_ID}"
