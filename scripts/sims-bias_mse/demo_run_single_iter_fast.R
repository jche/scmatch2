#!/usr/bin/env Rscript
# Demo: run a single iteration with fast methods and print a tiny summary
# Usage examples (from project root):
#   Rscript scripts/sims-bias_mse/demo_run_single_iter_fast.R
#   Rscript scripts/sims-bias_mse/demo_run_single_iter_fast.R hainmueller 1 "diff,bal1,ps_lm"
#   Rscript scripts/sims-bias_mse/demo_run_single_iter_fast.R toy 3 diff bal1 or_lm ps_lm

args <- commandArgs(trailingOnly = TRUE)

# Defaults: 1 DGP, 1 iteration, fast/stable methods
sim <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "toy"
id  <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 1L

# parse optional methods (space- or comma-separated). if none given, use fast set
parse_methods <- function(x) {
  if (length(x) == 0) return(character(0))
  m <- unlist(strsplit(x, split = ","))
  m <- trimws(tolower(m))
  m[nzchar(m)]
}
methods <- if (length(args) >= 3) parse_methods(args[-(1:2)]) else c("diff","bal1","or_lm","ps_lm")

cat("\n=== DEMO: run_single_iteration ===\n")
cat(sprintf("sim        : %s\n", sim))
cat(sprintf("iteration  : %d\n", id))
cat(sprintf("methods    : %s\n", if (length(methods)) paste(methods, collapse = ", ") else "ALL"))
cat("==================================\n\n")

# Build the command; prefer system2 with vector args for safety
runner <- "scripts/sims-bias_mse/run_single_iteration.R"
if (!file.exists(runner)) {
  stop(sprintf("Cannot find %s. Run this demo from the project root.", runner))
}

cmd_args <- c(runner, sim, as.character(id), methods)
status <- system2("Rscript", args = cmd_args)

if (!identical(status, 0L)) {
  stop(sprintf("run_single_iteration.R did not exit cleanly (status %s).", status))
}

# Where the script writes results:
out_file <- file.path("data", "outputs", "sims-bias_mse", sim, sprintf("iter_%04d.csv", id))
if (!file.exists(out_file)) {
  stop(sprintf("Result file not found: %s\nCheck the output_dir in run_single_iteration.R.", out_file))
}

# Quick peek at the result row
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

df <- read_csv(out_file, show_col_types = FALSE)
if (nrow(df) != 1) warning("Expected a single row in the result CSV.")

# Show the requested methods + a few meta columns
show_cols <- unique(c("sim_type","runid","seed","true_ATT", methods))
show_cols <- intersect(show_cols, names(df))

cat("\n--- Result summary ---\n")
print(df %>% select(all_of(show_cols)))

# Sanity: at least one requested method should be finite
req <- intersect(methods, names(df))
any_finite <- length(req) > 0 && any(is.finite(unlist(df[req])))
cat(sprintf("\nAt least one finite estimate among requested methods: %s\n",
            if (any_finite) "YES" else "NO"))
cat(sprintf("Result file: %s\n", out_file))
cat("\nDone.\n")
