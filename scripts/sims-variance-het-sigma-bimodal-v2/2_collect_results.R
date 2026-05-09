#!/usr/bin/env Rscript
# scripts/sims-variance-het-sigma-bimodal-v2/2_collect_results.R
#
# Combines individual .rds files from sims-variance-het-sigma-bimodal-v2
# into a single combined_results.rds / combined_results.csv.
#
# Usage:
#   Rscript 2_collect_results.R [output_name]
#   output_name defaults to "sims-variance-het-sigma-bimodal-v2"

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
})

args        <- commandArgs(trailingOnly = TRUE)
output_name <- if (length(args) >= 1) args[[1]] else "sims-variance-het-sigma-bimodal-v2"

out_base <- here::here("data/outputs", output_name)
paths <- list(
  out_base       = out_base,
  individual_dir = file.path(out_base, "individual"),
  combined_csv   = file.path(out_base, "combined_results.csv"),
  combined_rds   = file.path(out_base, "combined_results.rds")
)

dir.create(paths$out_base,       showWarnings = FALSE, recursive = TRUE)
dir.create(paths$individual_dir, showWarnings = FALSE, recursive = TRUE)

cat("Collecting .rds files from:\n  ", paths$individual_dir, "\n")

files <- list.files(
  path       = paths$individual_dir,
  pattern    = "^bimodal_v2_iter_[0-9]+_overlap_.*_s1extra_.*[.]rds$",
  full.names = TRUE
)

if (length(files) == 0) stop("No files found in: ", paths$individual_dir)
cat("Found ", length(files), " files.\n", sep = "")

read_one <- function(fp) {
  out <- tryCatch(readRDS(fp), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  out <- as_tibble(out)
  out$source_file <- basename(fp)
  out
}

res_list <- lapply(files, read_one)
failed   <- files[vapply(res_list, is.null, logical(1))]
combined <- bind_rows(Filter(Negate(is.null), res_list))

cat("Combined rows: ", nrow(combined), "\n", sep = "")
if (length(failed) > 0) {
  cat("Failed to read ", length(failed), " files:\n", sep = "")
  cat(paste("  -", basename(failed)), sep = "\n")
  cat("\n")
}

saveRDS(combined, paths$combined_rds)
readr::write_csv(combined, paths$combined_csv)

cat("Saved:\n")
cat("  RDS: ", paths$combined_rds, "\n", sep = "")
cat("  CSV: ", paths$combined_csv, "\n", sep = "")
