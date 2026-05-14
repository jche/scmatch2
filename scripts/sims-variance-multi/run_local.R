# scripts/sims-variance-multi/run_local.R
#
# Run multiple iterations locally in parallel using furrr.
# Useful for smoke-testing and quick sanity checks before submitting
# the full array job to SLURM.
#
# Usage (interactive):
#   source(here::here("scripts/sims-variance-multi/run_local.R"))
#
# Usage (command line):
#   Rscript scripts/sims-variance-multi/run_local.R [n_iter] [output_name]
#
#   n_iter       Number of iterations to run (default: 10)
#   output_name  Subdirectory under data/outputs/ (default: sims-variance-multi)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(furrr)
  library(tictoc)
  library(here)
})

# ── Configuration ─────────────────────────────────────────────────────────────
.args       <- commandArgs(trailingOnly = TRUE)
N_ITER      <- if (length(.args) >= 1) as.integer(.args[[1]]) else 100L
OUTPUT_NAME <- if (length(.args) >= 2) .args[[2]] else "sims-variance-multi"
SAVE_OUTPUT <- TRUE

# ── Load utilities ────────────────────────────────────────────────────────────
source(here::here("scripts/sims-variance-multi/0_utils.R"))

if (SAVE_OUTPUT) {
  paths <- get_sim_paths(OUTPUT_NAME)
  message(sprintf("Saving output to: %s", paths$individual_dir))
}

# ── Optionally save ───────────────────────────────────────────────────────────
save_path = NULL
if (SAVE_OUTPUT) {
  save_path = paths$individual_dir
}



# ── Run in parallel ───────────────────────────────────────────────────────────
n_workers <- parallel::detectCores() - 2
message(sprintf(
  "Running %d iteration(s) × %d design cells in parallel (%d workers)  [save=%s]",
  N_ITER, nrow(DESIGN_GRID), n_workers, SAVE_OUTPUT
))

plan(multisession, workers = n_workers)
tic()

all_res <- future_map( 1:N_ITER,
                       .f = sim_master_multi,
                       save_path = save_path,
                       overwrite = FALSE,
                       .options = furrr_options(seed = NULL),
                       .progress = TRUE
)

toc()
plan(sequential)



cat( "Simulation complete\n" )


