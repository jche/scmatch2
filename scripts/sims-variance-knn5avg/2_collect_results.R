# scripts/sims-variance/2_collect_results.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(here)
})

# Where your toy_match_infer saves individual results
paths <- list(
  out_base       = here::here("data/outputs/sims-variance-knn5avg"),
  individual_dir = here::here("data/outputs/sims-variance-knn5avg/individual"),
  combined_csv   = here::here("data/outputs/sims-variance-knn5avg/combined_results.csv"),
  combined_rds   = here::here("data/outputs/sims-variance-knn5avg/combined_results.rds")
)


dir.create(paths$out_base, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$individual_dir, showWarnings = FALSE, recursive = TRUE)

cat("Collecting .rds files from:\n  ", paths$individual_dir, "\n")

files <- list.files(
  path = paths$individual_dir,
  pattern = "^toy_match_infer_iter_\\d+_overlap_.*_error_.*\\.rds$",
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
failed <- files[vapply(res_list, is.null, logical(1))]

combined <- bind_rows(res_list)

cat("Combined rows: ", nrow(combined), "\n", sep = "")
if (length(failed) > 0) {
  cat("Failed to read ", length(failed), " files:\n", sep = "")
  cat(paste("  -", basename(failed)), sep = "\n")
  cat("\n")
}

# Save outputs
saveRDS(combined, paths$combined_rds)
readr::write_csv(combined, paths$combined_csv)

cat("Saved:\n")
cat("  RDS: ", paths$combined_rds, "\n", sep = "")
cat("  CSV: ", paths$combined_csv, "\n", sep = "")
