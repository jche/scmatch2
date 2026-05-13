# scripts/sims-variance-multi/2_collect.R
#
# Collect per-iteration .rds files into a single combined CSV + RDS.
#
# Usage:
#   Rscript 2_collect.R [output_name]
#   output_name  default: sims-variance-multi

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(here)
})

args        <- commandArgs(trailingOnly = TRUE)
output_name <- if (length(args) >= 1) args[[1]] else "sims-variance-multi"

source(here::here("scripts/sims-variance-multi/0_utils.R"))
paths <- get_sim_paths(output_name)

cat("Collecting .rds files from:\n  ", paths$individual_dir, "\n")

files <- list.files(
  path       = paths$individual_dir,
  pattern    = "^iter_\\d+\\.rds$",
  full.names = TRUE
)

if (length(files) == 0)
  stop("No files found in: ", paths$individual_dir)

cat("Found ", length(files), " files.\n", sep = "")

read_one <- function(fp) {
  out <- tryCatch(readRDS(fp), error = function(e) NULL)
  if (is.null(out)) { warning("Could not read: ", fp); return(NULL) }
  as_tibble(out)
}

res_list <- map(files, read_one)
failed   <- files[map_lgl(res_list, is.null)]
combined <- bind_rows(res_list)

cat("Combined rows: ", nrow(combined), "\n", sep = "")
cat("Replications:  ", n_distinct(combined$runID), "\n", sep = "")
cat("Design cells:  ", n_distinct(combined$overlap_label,
                                   combined$error_type,
                                   combined$sigma1_extra), "\n", sep = "")

if (length(failed) > 0) {
  cat("Failed to read ", length(failed), " files:\n", sep = "")
  cat(paste("  -", basename(failed), collapse = "\n"), "\n")
}

saveRDS(combined, paths$combined_rds)
write_csv(combined,  paths$combined_csv)

cat("Saved:\n  RDS: ", paths$combined_rds,
    "\n  CSV: ", paths$combined_csv, "\n")
