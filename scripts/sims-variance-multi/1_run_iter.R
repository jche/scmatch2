# scripts/sims-variance-multi/1_run_iter.R
#
# Unified entrypoint — works as a SLURM array job, a plain Rscript call,
# or an interactive source() in RStudio/R console.
#
# ── SLURM usage ───────────────────────────────────────────────────────────────
#   sbatch run_slurm.sh        (sets SLURM_ARRAY_TASK_ID automatically)
#
# ── Command-line usage ────────────────────────────────────────────────────────
#   Rscript 1_run_iter.R <iter_id> [output_name]
#
# ── Interactive / local usage ─────────────────────────────────────────────────
#   source(here::here("scripts/sims-variance-multi/1_run_iter.R"))
#
#   When no SLURM env var and no command-line args are present the script
#   falls back to local defaults (ITER = 1, no file save) and prints
#   diagnostic summaries at the end.  Adjust ITER / OUTPUT_NAME / SAVE_OUTPUT
#   in the LOCAL DEFAULTS block below as needed.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

# ── Mode detection ────────────────────────────────────────────────────────────
.slurm_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
.args     <- commandArgs(trailingOnly = TRUE)

verbose = FALSE
if (nchar(.slurm_id) > 0) {
  # ── SLURM array job ─────────────────────────────────────────────────────────
  iter        <- as.integer(.slurm_id)
  output_name <- if (length(.args) >= 1) .args[[1]] else "sims-variance-multi"
  save_output <- TRUE
  is_local    <- FALSE
  .mode       <- "SLURM"

} else if (length(.args) >= 1) {
  # ── Plain Rscript call ──────────────────────────────────────────────────────
  iter        <- as.integer(.args[[1]])
  output_name <- if (length(.args) >= 2) .args[[2]] else "sims-variance-multi"
  save_output <- TRUE
  is_local    <- FALSE
  .mode       <- "Rscript"

} else {
  cat( "Local mode!" )
  # ── LOCAL DEFAULTS (interactive / source()) ─────────────────────────────────
  iter        <- 1L           # which replication to run
  output_name <- "sims-variance-multi"
  save_output <- FALSE        # set TRUE to write iter_NNNN.rds to data/outputs/
  is_local    <- TRUE
  .mode       <- "local"
  verbose = TRUE
}

if (is.na(iter) || iter <= 0L)
  stop("iter_id must be a positive integer", call. = FALSE)

# ── Load utilities ────────────────────────────────────────────────────────────
source(here::here("scripts/sims-variance-multi/0_utils.R"))

# ── Run ───────────────────────────────────────────────────────────────────────
message(sprintf("[%s] iter=%d  output=%s  cells=%d",
                .mode, iter, output_name, nrow(DESIGN_GRID)))
t0  <- Sys.time()

res <- sim_master_multi(iteration = iter, verbose=verbose)

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
message(sprintf("Done in %.1fs  (%d rows)", elapsed, nrow(res)))

# ── Optionally save ───────────────────────────────────────────────────────────
if (save_output) {
  paths    <- get_sim_paths(output_name)
  out_file <- file.path(paths$individual_dir, sprintf("iter_%04d.rds", iter))
  saveRDS(res, out_file)
  message("Saved: ", out_file,
          "  (", nrow(res), " rows, ",
          round(sum(res$time_secs, na.rm = TRUE), 1), "s total)")
}



# ── Local diagnostics (only when run interactively) ───────────────────────────
if (is_local) {
  cat("\n── Design cells run ──\n")
  res %>%
    distinct(overlap_label, error_type, sigma1_extra, k_match) %>%
    arrange(overlap_label, error_type, sigma1_extra, k_match) %>%
    print(n = Inf)

  cat("\n── Status counts ──\n")
  print(table(res$status))

  cat("\n── Coverage by estimator × error_type × common_label × k_match ──\n")
  res %>%
    filter(!is.na(CI_lower), !is.na(SATT)) %>%
    mutate(covered = CI_lower <= SATT & SATT <= CI_upper) %>%
    group_by(inference_method, error_type, common_label, k_match) %>%
    summarise(
      n        = n(),
      coverage = mean(covered),
      mean_SE  = mean(SE, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    print(n = Inf)

  cat("\n── ATT estimates (first few rows) ──\n")
  print(head(res %>% select(runID, inference_method, overlap_label, error_type,
                             common_label, k_match, att_est, SE, SATT), 20))
}
