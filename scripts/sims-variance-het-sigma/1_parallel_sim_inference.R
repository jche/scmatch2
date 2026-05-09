# scripts/sims-variance-het-sigma/1_parallel_sim_inference.R
#
# SLURM usage:
#   sbatch run-adaptive-k2-het-sigma.sh
#
# Direct usage:
#   Rscript 1_parallel_sim_inference.R <iter_id> [output_name]
#
# Each job loops over all 5 overlap levels × 2 sigma1_extra values.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 1_parallel_sim_inference.R <iter_id> [output_name]",
       call. = FALSE)
}

iter        <- as.integer(args[[1]])
output_name <- if (length(args) >= 2) args[[2]] else "sims-variance-het-sigma"

if (is.na(iter) || iter <= 0) stop("iter_id must be a positive integer", call. = FALSE)

source(here::here("scripts/sims-variance-het-sigma/0_sim_inference_utils.R"))

paths <- get_sim_paths(output_name)
dir.create(paths$individual_dir, showWarnings = FALSE, recursive = TRUE)

overlap_labels <- c("very_low", "low", "mid", "high", "very_high")
sigma1_extras  <- c(0, 0.5)

message("Starting iter=", iter,
        " output=", output_name,
        " overlaps={", paste(overlap_labels, collapse = ","), "}",
        " sigma1_extras={", paste(sigma1_extras, collapse = ","), "}")

for (overlap_label in overlap_labels) {
  for (sigma1_extra in sigma1_extras) {

    res <- sim_master_het(
      iteration    = iter,
      overlap_label = overlap_label,
      sigma1_extra = sigma1_extra,
      rad_method   = "adaptive",
      k_match      = 2,
      est_method   = "scm",
      K_tt         = 2
    )

    sigma_tag <- sprintf("s1extra_%.2f", sigma1_extra)
    out_file  <- file.path(
      paths$individual_dir,
      sprintf("het_sigma_iter_%04d_overlap_%s_%s.rds", iter, overlap_label, sigma_tag)
    )

    saveRDS(res, out_file)
    message("Saved: ", out_file)
  }
}

message("Done iter=", iter)
