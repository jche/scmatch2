# Scripts Directory

This directory contains all scripts used to generate the figures, simulations, and empirical analyses for:

> **Caliper Synthetic Matching: Generalized Radius Matching with Local Synthetic Controls**
> Che, Meng, and Miratrix.  (2025).

The paper introduces CSM, a matching method that combines adaptive calipers with local synthetic control reweighting. The scripts here reproduce every figure and result in the paper and its online supplement.

---

## Directory Overview

```
scripts/
├── lib/                       # Shared function libraries (sourced by other scripts)
├── figs/                      # Generated figures
├── sims-bias_mse/             # Bias/MSE simulations (Section 5)
├── sims-variance-multi/       # Multifactor variance simulation (current paper, 5×2×2 design)
├── sims-variance/             # Variance simulations (earlier, overlap-only; superseded)
├── sims-variance-knn5avg/     # Variant: knn5avg aggregation (robustness check)
├── ferman-analysis/           # Empirical application: Brazilian education data (Section 6)
├── lalonde-analysis/          # Empirical application: LaLonde data (Online Supplement)
├── datagen/                   # Data generation utilities
├── demo/                      # Package demos and usage examples
├── boot/                      # Bootstrap inference experiments (exploratory)
└── draft-inference-scripts/   # Earlier draft inference scripts (superseded)
```

---

## Subfolders

### `figs/` — Paper figures

Scripts that produce the figures appearing in the paper and supplement. Run these after the simulations and analyses have been completed.

| Script | Figure | Paper location |
|---|---|---|
| `fig_toy_example.R` | Toy examples (joint balance motivation) | Section 2, Figure 1 |
| `fig_overlap_clusters_1d.R` | 1D overlap/caliper illustration | Section 3 |
| `fig_cem_vs_csm.R` | CEM vs. CSM comparison | Section 3 |
| `fig-sim_toy_3_overlaps.R` | Toy simulation across overlap levels | Section 5 |
| `fig_sim_toy.R` | Toy simulation summary figure | Section 5 |
| `fig_sim_canonicals.R` | Canonical DGP results (Kang, Hainmueller, ACIC) | Section 5 |
| `fig_clt-verification.R` | CLT verification for variance estimator | Online Supplement |
| `toy_demo.R` | Additional toy demo figure | Online Supplement |

Most figure scripts read pre-computed outputs from `data/outputs/` and write PNGs to `figures/`. `generate_toy_background_df.R` is a small helper that pre-computes the background density grid for the toy figures; its output is cached as `data/inputs/toy_background_df.rds` and regenerated automatically if missing.

---

### `sims-bias_mse/` — Bias and MSE simulations (Section 5)

The main simulation study comparing CSM to existing matching and weighting methods on four data-generating processes (DGPs): `toy`, `kang`, `hainmueller`, and `acic`.

See [`sims-bias_mse/README_simulations.md`](sims-bias_mse/README_simulations.md) for full workflow details. The typical pipeline is:

1. **Quick demo:** `Rscript sims-bias_mse/demo_run_single_iter_fast.R`
2. **Full run (local):** `Rscript sims-bias_mse/run_local_all.R <n_iterations>`
3. **Full run (HPC/SLURM):** `sbatch sims-bias_mse/run_slurm_all.sh`
4. **Collect results:** `Rscript sims-bias_mse/collect_sim_results.R`
5. **Plot:** run `figs/fig_sim_canonicals.R` and `figs/fig_sim_toy.R`

Key scripts:

- `run_single_iteration.R` — core worker script; runs one DGP iteration
- `run_local_all.R` — runs all DGPs sequentially on a local machine
- `run_slurm_all.sh` — SLURM array job for HPC
- `collect_sim_results.R` — combines per-iteration CSVs into a single file

---

### `sims-variance/` — Variance and inference simulations

Simulates coverage and performance of the variance estimator from Meng (2025) applied to CSM, across multiple overlap regimes and error structures.

See [`sims-variance/sims-variance-replication.md`](sims-variance/sims-variance-replication.md) for full workflow details. The numbered scripts form the pipeline:

1. `0_sim_inference_utils.R` — shared utilities (`sim_master()`, `toy_match_infer()`)
2. `1_parallel_sim_inference.R` — SLURM task entrypoint (one iteration)
3. `2_collect_results.R` — collects `.rds` outputs into a single CSV
4. `3_summarize_results.R` — produces summary tables
5. `4_debug_one_iter_V_vs_VE.R` — diagnostic for a single iteration

Run script: `run-variance-sims.sh`

---

### `sims-variance-multi/` — **Multifactor variance simulation (current paper)**

The definitive 3-factor simulation evaluating all four variance estimators. Supersedes `sims-variance/` and `sims-variance-het-sigma-bimodal-v2/`.

**Factors** (5 × 2 × 2 = 20 cells, 998 replications):

| Factor | Levels |
|---|---|
| Overlap | very_low, low, mid, high, very_high |
| Error structure | homo (σ₀ = 0.5 constant), het (bimodal σ₀) |
| Common variance | common (σ₁ = σ₀), no_common (σ₁ = σ₀ + 2) |

**Estimators:** `homo`, `het`, `alt_common`, `alt_tt`

**Performance measure:** coverage of τ_SATT at 95%.

Pipeline:

1. `0_utils.R` — DGP, estimators, `one_iter()`, `sim_master_multi()`
2. `1_run_iter.R` — SLURM entrypoint: one job runs one iteration × all 20 cells
3. `2_collect.R` — collects per-iteration `.rds` files into a combined CSV
4. `3_summarize.R` — computes coverage (with MCSEs via `simhelpers`) and saves the faceted coverage plot

Run script: `run_slurm.sh` (array 1–998)

---

### `sims-variance/` — Variance simulations (earlier, overlap-only)

Earlier simulation varying only overlap (homoskedastic DGP, two estimators). Superseded by `sims-variance-multi/` but retained for reference.

---

### `ferman-analysis/` — Empirical application: Brazilian education data (Section 6)

Replicates the within-study comparison using the Ferman (2021) data on Brazil's "Jovem de Futuro" program. Data are in `ferman-analysis/data/`.

Pipeline (numbered scripts):

1. `01-clean-data-for-analysis.R` — loads and cleans `Final.dta`; saves `ferman_for_analysis.rds`
2. `02-core-ferman-analysis.R` — runs CSM matching and estimates ATT
3. `03-ferman-figures.R` — produces figures for Section 6

Additional figure scripts (each generates one paper figure):

- `fig-hist-top-k-distances.R` — histogram of top-k match distances
- `fig-love-plot-ferman.R` — Love plot of covariate balance
- `fig-table-ESS.R` — effective sample size table
- `marginal-mean-table.R` — marginal means by treatment group
- `n-feasible-table.R` — number of feasible matched units

The original Ferman (2021) Stata code is preserved in `2021-Ferman-Stata-Code/` for reference.

---

### `lalonde-analysis/` — Empirical application: LaLonde data (Online Supplement)

Applies CSM to the classic LaLonde (1986) dataset using the CPS comparison group. Data are in `lalonde-analysis/data/` and `data/inputs/`.

Pipeline:

1. `01-clean-data-lalonde.R` — loads `lalonde_w_cps.RData`, sets up treatment/control coding
2. `02-core-lalonde-analysis.R` — runs CSM and computes ATT estimates
3. `03-lalonde-figures.R` — produces figures for the online supplement

Older exploratory scripts (`old_lalonde_code.R`, `twang_exploration.R`, etc.) are retained but not part of the main pipeline.

---

### `lib/` — Shared function libraries

Scripts that define functions sourced by multiple other scripts. Nothing here is meant to be run directly.

| File | Purpose | Sourced by |
|---|---|---|
| `wrappers.R` | `get_att_*()` estimator wrappers; `parse_form()` | sims, figs, tests |
| `sim_runner.R` | `run_all_methods()`, `expand_method_groups()`, `ALL_METHODS` | sims, tests |
| `plot_sim.R` | Plotting functions for simulation results (`acic_plot_sim_type()` etc.) | `figs/fig_sim_*.R` |
| `plot_toy.R` | Plotting functions for toy examples | `figs/fig_toy_example.R` |

---

### `datagen/` — Data generation utilities

Scripts that construct input datasets used by simulations and demos:

- `gen_lalonde_w_cps.R` — constructs the LaLonde + CPS dataset
- `gen_six_points.R` — generates the six-point toy dataset
- `gen_uniform.R` — generates a uniform covariate dataset

---

### `demo/` — Package usage examples

Standalone scripts illustrating the `SCMatch` package API. Useful for learning the package or verifying an installation.

- `demo.R` — basic CSM matching demo
- `example-lalonde.R` — end-to-end LaLonde example
- `CSM-vs-CEM.R` — side-by-side comparison of CSM and CEM
- `s3-learning.R` — exploration of the package's S3 class structure

---

### `boot/` — Bootstrap inference experiments (exploratory)

Scripts exploring bootstrap-based confidence intervals as an alternative to the analytic variance estimator. These are not part of the final paper but informed the decision to use the Meng (2025) variance estimator instead.

- `main.R` — top-level driver
- `block-bootstrap.R`, `spatial-block-sampling.R`, `subsampling.R` — bootstrap variants
- `boot_CSM_simulation_code.R` — simulation comparing bootstrap CI methods
- `bootstrap_simulation_runs.R` — runs bootstrap simulations
- `figures/` — CI comparison plots
- `output/` — saved simulation results
- `development-otsu/` — development work on the Otsu-Rai bootstrap method

---

### `draft-inference-scripts/` — Earlier inference scripts (superseded)

An earlier version of the inference simulation framework, written during initial development. Superseded by `sims-variance/`. Retained for reference.

Key files: `ferman-experiments.R`, `utils-replicate-ferman.R`, `verify_CLT.R`.

---

## Data flow

```
datagen/          →  data/inputs/
sims-bias_mse/    →  data/outputs/sims-bias_mse/   →  figs/fig_sim_*.R  →  figures/
sims-variance/    →  data/outputs/sims-variance/   →  figs/fig_clt-*.R  →  figures/
ferman-analysis/  →  (in-memory)                   →  figs/fig-*ferman* →  figures/
lalonde-analysis/ →  (in-memory)                   →  figs/ (OS)        →  figures/
```

All scripts use `here::here()` for paths and assume the project root is the working directory (open `scmatch2.Rproj` in RStudio, or set working directory to the repo root before running).
