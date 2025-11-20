#!/bin/bash
#SBATCH -J all_sims
#SBATCH -o /homes2/xmeng/scmatch2/logs/all_%A_%a.out
#SBATCH -e /homes2/xmeng/scmatch2/logs/all_%A_%a.err
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 02:30:00
#SBATCH --array=1-2

# ---------------------------------------------------------------------
# SLURM job array: run all four DGPs (acic, hainmueller, kang, toy)
# Output written to data/outputs/sims/<sim>/iter_####.csv
# ---------------------------------------------------------------------

project_dir="/homes2/xmeng/scmatch2/"
module load R/4.3.2
mkdir -p "${project_dir}/logs"
cd "${project_dir}" || exit 1

ITER=${SLURM_ARRAY_TASK_ID}

echo "========================================="
echo "Job ${SLURM_JOB_ID}  | Iter ${ITER}"
echo "Start: $(date)"
echo "Methods: ${METHODS:-ALL}"
echo "========================================="

IFS=', ' read -r -a METHODS_ARR <<< "${METHODS}"

run_one() {
  local sim="$1"
  echo "--- ${sim} / iter ${ITER} ---"
  if [[ -z "${METHODS}" ]]; then
    Rscript scripts/sims/run_single_iteration.R "${sim}" "${ITER}"
  else
    Rscript scripts/sims/run_single_iteration.R "${sim}" "${ITER}" "${METHODS_ARR[@]}"
  fi
}

set -e
for sim in acic hainmueller kang toy; do
  run_one "${sim}"
done
set +e

echo "Done: $(date)"
