# scripts/sims-variance-multi/run_local.R
#
# Run multiple iterations locally (no SLURM).  Useful for smoke-testing
# and quick sanity checks before submitting the full array job.
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
  library(here)
})

# ── Configuration ─────────────────────────────────────────────────────────────
.args       <- commandArgs(trailingOnly = TRUE)
N_ITER      <- if (length(.args) >= 1) as.integer(.args[[1]]) else 10L
OUTPUT_NAME <- if (length(.args) >= 2) .args[[2]] else "sims-variance-multi"
SAVE_OUTPUT <- TRUE #length(.args) >= 3 && .args[[3]] == "save"

# ── Load utilities ────────────────────────────────────────────────────────────
source(here::here("scripts/sims-variance-multi/0_utils.R"))

if (SAVE_OUTPUT) {
  paths <- get_sim_paths(OUTPUT_NAME)
  message(sprintf("Saving output to: %s", paths$individual_dir))
}

# ── Run ───────────────────────────────────────────────────────────────────────
message(sprintf(
  "Running %d iteration(s) × %d design cells  [save=%s]",
  N_ITER, nrow(DESIGN_GRID), SAVE_OUTPUT
))

t_total <- Sys.time()
all_res  <- vector("list", N_ITER)

for (i in seq_len(N_ITER)) {
  t0  <- Sys.time()
  res <- sim_master_multi(iteration = i)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message(sprintf("  iter %3d / %d  done in %.1fs  (%d rows, %d errors)",
                  i, N_ITER, elapsed, nrow(res),
                  sum(res$status != "ok", na.rm = TRUE)))

  if (SAVE_OUTPUT) {
    out_file <- file.path(paths$individual_dir, sprintf("iter_%04d.rds", i))
    saveRDS(res, out_file)
  }

  all_res[[i]] <- res
}

elapsed_total <- round(as.numeric(difftime(Sys.time(), t_total, units = "secs")), 1)
message(sprintf("\nAll %d iterations done in %.1fs (%.1fs / iter)",
                N_ITER, elapsed_total, elapsed_total / N_ITER))

# ── Combine ───────────────────────────────────────────────────────────────────
results <- bind_rows(all_res)
