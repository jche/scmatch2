# scripts/sims-variance-multi/2_collect.R
#
# Collect per-iteration .rds files into a single combined CSV + RDS.
# Logic lives in collect_sim_results() in 0_utils.R.
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

collect_sim_results(output_name)
