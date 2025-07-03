for i in {1..500}; do
  if [ ! -f "/homes2/xmeng/scmatch2/data/outputs/inf_toy/individual/results_iter_${i}.rds" ]; then
    echo "$i"
  fi
done > missing_ids.txt

