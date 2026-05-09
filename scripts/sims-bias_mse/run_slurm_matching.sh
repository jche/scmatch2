#!/bin/bash
#SBATCH -J matching_sims
#SBATCH -o /homes2/xmeng/scmatch2/logs/matching_%A_%a.out
#SBATCH -e /homes2/xmeng/scmatch2/logs/matching_%A_%a.err
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --array=1-1000

# ---------------------------------------------------------------------
# SLURM job array: run CSM + CEM methods only (full-bin and half-bin)
# for all four DGPs (acic, hainmueller, kang, toy).
#
# Methods run: matching group  = csm_scm, csm_avg, cem_scm, cem_avg, onenn
#              matching_half   = csm_scm_half, csm_avg_half, cem_scm_half, cem_avg_half
#
# Submit:
#   sbatch scripts/sims-bias_mse/run_slurm_matching.sh
#
# To run a subset of sims, edit the loop below or pass --export:
#   sbatch --export=ALL,SIMS="toy" scripts/sims-bias_mse/run_slurm_matching.sh
# ---------------------------------------------------------------------

project_dir="/homes2/xmeng/scmatch2"
module load R/4.3.2
mkdir -p "${project_dir}/logs"
cd "${project_dir}" || exit 1

ITER=${SLURM_ARRAY_TASK_ID}
METHODS="matching matching_half"

echo "========================================="
echo "Job ${SLURM_JOB_ID}  | Iter ${ITER}"
echo "Start: $(date)"
echo "Methods: ${METHODS}"
echo "========================================="

run_one() {
  local sim="$1"
  echo "--- ${sim} / iter ${ITER} ---"
  Rscript scripts/sims-bias_mse/run_single_iteration.R "${sim}" "${ITER}" ${METHODS}
}

set -e
for sim in ${SIMS:-acic hainmueller kang toy}; do
  run_one "${sim}"
done
set +e

echo "Done: $(date)"
