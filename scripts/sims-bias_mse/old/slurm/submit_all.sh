#!/bin/bash
# Submit all three simulations

project_dir="/homes2/xmeng/scmatch2/"

cd ${project_dir}

echo "Submitting ACIC simulation..."
sbatch scripts/slurm/submit_acic.sh

echo "Submitting Hainmueller simulation..."
sbatch scripts/slurm/submit_hainmueller.sh

echo "Submitting Kang simulation..."
sbatch scripts/slurm/submit_kang.sh

echo ""
echo "All jobs submitted!"
echo "Check status with: squeue -u $USER"
