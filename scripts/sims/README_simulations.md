# Simulation Workflow Guide

This folder (`scripts/sims/`) contains all scripts used to reproduce the simulation experiments in **CSM**.

The workflow supports three usage modes:

1. **Fast single-run demo** – quick smoke test on one DGP.  
2. **Local multi-run** – run all DGPs sequentially on your workstation.  
3. **Cluster (SLURM) run** – run all DGPs in parallel on an HPC cluster.

---

## ⚙️ Folder structure

```
scripts/sims/
├── run_single_iteration.R      # core script: runs one iteration for one DGP
├── run_local_all.R             # runs all DGPs sequentially (local CPU)
├── run_slurm_all.sh            # SLURM array job (for HPC)
├── demo_run_single_iter_fast.R # quick demo for fast methods (local)
└── ...
```

Output directories:

```
data/outputs/sims/<sim_type>/
  ├── iter_0001.csv             # per-iteration results
  ├── iter_0002.csv
  └── ...
```

where `<sim_type>` ∈ `{acic, hainmueller, kang, toy}`.

---

## 🚀 1. Fast demo run

For a minimal sanity check (≈ 10 seconds):

```
Rscript scripts/sims/demo_run_single_iter_fast.R
```

This runs the `toy` DGP for one iteration using only fast estimators  
(`diff`, `bal1`, `or_lm`, `ps_lm`) and prints a one-line summary.

Alternatively, **you can open the script in RStudio and run all lines interactively**  
to see the workflow step by step.

Output → `data/outputs/sims/toy/iter_0001.csv`.

---

## 🧩 2. Local multi-run (no cluster)

Run several iterations for all DGPs sequentially:

```
Rscript scripts/sims/run_local_all.R 5
```

- `5` = number of iterations (default = 5).  
- Optional list of methods (space or comma separated):

```
Rscript scripts/sims/run_local_all.R 3 diff bal1 or_lm ps_lm
```

Results for each DGP are saved to  
`data/outputs/sims/<sim_type>/iter_####.csv`.

---

## 🧮 3. Full SLURM run (cluster)

Submit the array job:

```
cd scripts/sims
sbatch run_slurm_all.sh
```

Each array element runs all four DGPs once (e.g., iterations 1–1000).

### Logs

Logs are written to:

```
/homes2/<your_username>/scmatch2/logs/all_<jobid>_<arrayid>.out
/homes2/<your_username>/scmatch2/logs/all_<jobid>_<arrayid>.err
```

To monitor jobs:

```
squeue -u $USER
```

To view a log:

```
less /homes2/<your_username>/scmatch2/logs/all_<jobid>_1.out
```

---

## 📦 4. Collect results (after SLURM run)

Combine iteration CSVs and print performance summaries:

```
Rscript scripts/analysis/collect_slurm_results.R
```

Alternatively, **you can open and run `scripts/analysis/collect_slurm_results.R`  
directly in RStudio** to interactively inspect summaries.

Creates combined files:

```
data/outputs/sims/<sim_type>_combined.csv
```

and copies summarized results to:

```
data/outputs/sim_canonical_results/
```

---

## 🖼️ 5. Generate figures

Once results are collected, you can reproduce the paper-quality plots via:

```
scripts/figs/fig_sim_canonicals.R   # canonical DGPs (ACIC, Hainmueller, Kang)
scripts/figs/fig_sim_toy.R          # toy simulation
```

These scripts read the combined CSVs and produce all simulation figures.

---

## ⚠️ 6. Important: Adjust your paths

Before running on your own system:

1. **Edit path variables** inside scripts:  
   - In SLURM script (`run_slurm_all.sh`), update:
     ```
     project_dir="/homes2/xmeng/scmatch2/"
     ```
     → change to your own absolute project path.

   - In the SBATCH header, confirm log paths:
     ```
     #SBATCH -o /homes2/xmeng/scmatch2/logs/all_%A_%a.out
     #SBATCH -e /homes2/xmeng/scmatch2/logs/all_%A_%a.err
     ```

2. Make sure the following directories exist (the scripts will create them if not):

```
data/outputs/sims/
logs/
```

3. Confirm that R and required packages load correctly:

```
module load R/4.3.2
Rscript -e "library(CSM); library(twang); library(kbal)"
```

---

## 🧠 7. Typical workflow summary

| Step | Command | Purpose |
|------|----------|----------|
| Quick check | `Rscript scripts/sims/demo_run_single_iter_fast.R` or run in RStudio | Verify environment |
| Local dev run | `Rscript scripts/sims/run_local_all.R 3` | Run small replication locally |
| Cluster run | `sbatch scripts/sims/run_slurm_all.sh` | Full-scale SLURM job array |
| Combine results | `Rscript scripts/analysis/collect_slurm_results.R` or run interactively | Summarize all outputs |
| Make figures | run `scripts/figs/fig_sim_canonicals.R` / `scripts/figs/fig_sim_toy.R` | Generate simulation plots |

---

## 💡 Tips

- You can limit methods on SLURM by setting an environment variable before submitting:

```
export METHODS="diff,bal1,or_lm,ps_lm"
sbatch scripts/sims/run_slurm_all.sh
```

- To debug interactively, try running one iteration manually:

```
Rscript scripts/sims/run_single_iteration.R hainmueller 1 diff bal1
```

---

**Author:** Xiang Meng  
**Last updated:** 2025-11-21  
**Contact:** xmeng@g.harvard.edu
