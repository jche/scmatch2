# scripts/sims-variance/1_parallel_sim_inference.R
#
# SLURM usage (via one of the three sh scripts):
#   sbatch run-adaptive-k5-scm.sh        # adaptive caliper, k=5, SCM weights
#   sbatch run-adaptive-k1-scm.sh        # adaptive caliper, k=1, SCM weights
#   sbatch run-knn-k5-avg.sh             # k-NN, k=5, average weights
#
# Direct usage:
#   Rscript 1_parallel_sim_inference.R <iter_id> <output_name> <rad_method> <k> <est_method>
#
# Arguments (positional):
#   iter_id      positive integer, SLURM_ARRAY_TASK_ID
#   output_name  subdirectory under data/outputs/  (e.g. sims-variance)
#   rad_method   "adaptive" or "knn"
#   k            integer caliper / neighbour count
#   est_method   "scm" or "average"

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

# ---- Parse args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript 1_parallel_sim_inference.R <iter_id> <output_name> <rad_method> <k> <est_method>",
    call. = FALSE
  )
}

iter       <- as.integer(args[[1]])
output_name <- args[[2]]
rad_method  <- args[[3]]
k_match     <- as.integer(args[[4]])
est_method  <- args[[5]]

if (is.na(iter) || iter <= 0) stop("iter_id must be a positive integer", call. = FALSE)
if (!rad_method %in% c("adaptive", "knn")) stop("rad_method must be 'adaptive' or 'knn'", call. = FALSE)
if (is.na(k_match) || k_match <= 0) stop("k must be a positive integer", call. = FALSE)
if (!est_method %in% c("scm", "average")) stop("est_method must be 'scm' or 'average'", call. = FALSE)

# ---- Source utils ----
source(here::here("scripts/sims-variance/0_sim_inference_utils.R"))

# ---- Paths ----
paths <- get_sim_paths(output_name)
dir.create(paths$individual_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Simulation config ----
N_val <- 600
grid_id <- 1

overlap_labels <- c("very_low", "low", "mid", "high", "very_high")
error_labels   <- c("homoskedastic")

# ---- Run one iteration across scenarios ----
message("Running iter=", iter,
        " output=", output_name,
        " rad_method=", rad_method,
        " k=", k_match,
        " est_method=", est_method,
        " across overlaps={", paste(overlap_labels, collapse = ","), "}",
        " errors={", paste(error_labels, collapse = ","), "}")

for (overlap_label in overlap_labels) {
  for (error_label in error_labels) {

    res <- sim_master(
      iteration  = iter,
      N          = N_val,
      overlap_label = overlap_label,
      error_label   = error_label,
      rad_method    = rad_method,
      k_match       = k_match,
      est_method    = est_method,
      k_dim  = 2,
      grid_id = grid_id
    )

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
