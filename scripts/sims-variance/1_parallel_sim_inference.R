# scripts/sims-variance/1_parallel_sim_inference.R
#
# SLURM usage (recommended):
#   sbatch --array=1-500 your_slurm_script.sh
# Each array task runs ONE iteration index (iter) across all overlap scenarios,
# saving one .rds per (iter, overlap_label, error_label).
#
# After the array finishes, run a separate combine script (or call
# collect_results_to_csv()).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

# ---- Parse iteration from command line ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 1_parallel_sim_inference.R <iter_id>", call. = FALSE)
}
iter <- as.integer(args[[1]])
if (is.na(iter) || iter <= 0) stop("iter_id must be a positive integer", call. = FALSE)

# iter <- 1

# ---- Source utils (adjusted directory) ----
source(here::here("scripts/sims-variance/0_sim_inference_utils.R"))

# ---- Paths ----
paths <- get_sim_paths()
dir.create(paths$individual_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Simulation config ----
# Fixed N is metadata; sim_master enforces nc=500, nt=100 for this replication
N_val <- 600
grid_id <- 1

# Overlap scenarios: controls uniform proportion determined inside sim_master via overlap_label
overlap_labels <- c("very_low", "low", "mid", "high", "very_high")

# For strict replication of the setting you described: homoskedastic only
error_labels <- c("homoskedastic")

# ---- Run one iteration across scenarios ----
message("Running iter=", iter,
        " across overlaps={", paste(overlap_labels, collapse = ","), "}",
        " errors={", paste(error_labels, collapse = ","), "}")

for (overlap_label in overlap_labels) {
  for (error_label in error_labels) {
    # overlap_label <- "very_low"; error_label <- "homoskedastic"

    res <- sim_master(
      iteration = iter,
      N = N_val,
      overlap_label = overlap_label,
      error_label = error_label,
      k_dim = 2,
      grid_id = grid_id
    )

    # Save per-scenario result; safe for parallel SLURM array
    out_file <- file.path(
      paths$individual_dir,
      sprintf(
        "toy_match_infer_iter_%04d_overlap_%s_error_%s.rds",
        iter, overlap_label, error_label
      )
    )

    saveRDS(res, out_file)
    message("Saved: ", out_file)
  }
}

message("Done iter=", iter)
