# scripts/sims-variance-multi/1_run_iter.R
#
# SLURM entrypoint — one array job = one replication × all 20 design cells.
#
# SLURM usage:
#   sbatch run_slurm.sh
#
# Direct usage (single replication, e.g. for testing):
#   Rscript 1_run_iter.R <iter_id> [output_name]
#
# Arguments:
#   iter_id      Positive integer — SLURM_ARRAY_TASK_ID (1 .. 998).
#   output_name  Subdirectory under data/outputs/ (default: sims-variance-multi).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

# ── Parse arguments ──────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
iter        <- as.integer(args[[1]])
output_name <- if (length(args) >= 2) args[[2]] else "sims-variance-multi"

if (is.na(iter) || iter <= 0L)
  stop("iter_id must be a positive integer", call. = FALSE)

# ── Load utilities ───────────────────────────────────────────────────────────
source(here::here("scripts/sims-variance-multi/0_utils.R"))

# ── Output paths ─────────────────────────────────────────────────────────────
paths <- get_sim_paths(output_name)

# ── Run ──────────────────────────────────────────────────────────────────────
message(sprintf("iter=%d  output=%s  cells=%d", iter, output_name, nrow(DESIGN_GRID)))

res <- sim_master_multi(iteration = iter)

# ── Save ─────────────────────────────────────────────────────────────────────
out_file <- file.path(
  paths$individual_dir,
  sprintf("iter_%04d.rds", iter)
)
saveRDS(res, out_file)
message("Saved: ", out_file,
        "  (", nrow(res), " rows, ",
        round(sum(res$time_secs, na.rm = TRUE), 1), "s total)")
