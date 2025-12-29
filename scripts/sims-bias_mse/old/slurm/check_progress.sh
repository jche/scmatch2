#!/bin/bash
# Check progress of SLURM simulations

project_dir="/homes2/xmeng/scmatch2/"

cd ${project_dir}

echo "=== Checking Simulation Progress ==="
echo ""

for sim in acic hainmueller kang; do
  dir="data/outputs/sim/${sim}"
  if [ -d "$dir" ]; then
    count=$(ls ${dir}/iter_*.csv 2>/dev/null | wc -l)
    echo "${sim}: ${count}/100 iterations complete"
  else
    echo "${sim}: directory not found"
  fi
done

echo ""
echo "=== Current SLURM Jobs ==="
squeue -u $USER

echo ""
echo "=== Recent Log Errors (last 5 lines from each type) ==="
for sim in acic hain kang; do
  latest_err=$(ls -t logs/${sim}_*.err 2>/dev/null | head -1)
  if [ -f "$latest_err" ]; then
    echo "--- ${sim} (from ${latest_err}) ---"
    tail -5 "$latest_err"
  fi
done
