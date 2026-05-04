#!/bin/bash
#SBATCH -J csm_het_sigma
#SBATCH -o logs/csm_het_sigma_%A_%a.out
#SBATCH -e logs/csm_het_sigma_%A_%a.err
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH --array=1-1000

project_dir="/homes2/xmeng/scmatch2"
r_script_path="${project_dir}/scripts/sims-variance-het-sigma/1_parallel_sim_inference.R"

module load R/4.3.2

mkdir -p "${project_dir}/logs"

if [ ! -f "$r_script_path" ]; then
  echo "Error: R script not found at $r_script_path"
  exit 1
fi

echo "Starting iter=${SLURM_ARRAY_TASK_ID} [het-sigma, k=2, adaptive, scm]"
cd "${project_dir}"

Rscript "${r_script_path}" "${SLURM_ARRAY_TASK_ID}" "sims-variance-het-sigma"

echo "Finished iter=${SLURM_ARRAY_TASK_ID}"
