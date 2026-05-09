#!/usr/bin/env Rscript
# Run a single iteration of the comprehensive simulation (+ toy)
#
# Usage:
# Rscript scripts/sims-bias_mse/canonical_run_single_iteration.R <sim_type> <iteration_id> [methods...]
#
# sim_type: one of 'acic','hainmueller','kang','toy'
# iteration_id: integer
# methods: optional list and/or comma-separated names or group shorthands
#          (see expand_method_groups() in scripts/sim_runner.R)

devtools::load_all()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript canonical_run_single_iteration.R <sim_type> <iteration_id> [methods...]")
}
sim_type     <- args[1]
iteration_id <- as.integer(args[2])

# --- parse optional methods ---------------------------------------------------
raw_methods <- character(0)
if (length(args) > 2) {
  margs       <- args[-(1:2)]
  raw_methods <- unlist(strsplit(margs, split = ","))
  raw_methods <- trimws(tolower(raw_methods))
  raw_methods <- raw_methods[nzchar(raw_methods)]
}

suppressPackageStartupMessages({
  require(tidyverse)
  require(mvtnorm)
  require(optweight)
  require(dbarts)
  require(tmle)
  require(AIPW)
  require(grf)
  require(twang)
  require(kbal)
  require(tictoc)
  require(CSM)
})
source(here::here("scripts/lib/wrappers.R"))
source(here::here("scripts/lib/sim_runner.R"))
source(here::here("R/utils.R"))
source(here::here("R/sim_data.R"))

# Expand group shorthands; default to all methods if none specified
SELECTED_METHODS <- if (length(raw_methods) == 0) ALL_METHODS else expand_method_groups(raw_methods)

cat(sprintf("\n=== Running %s simulation, iteration %d ===\n", sim_type, iteration_id))
cat("Methods: ", paste(SELECTED_METHODS, collapse = ", "), "\n")

output_dir <- file.path("data", "outputs", "sims-bias_mse-no-filter-unmatched", sim_type)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- seeds --------------------------------------------------------------------
seed_start <- 1000
if (sim_type == "acic") {
  set.seed(seed_start + iteration_id)
  current_seed <- sample(1:100000, 1)
  set.seed(current_seed)
} else {
  current_seed <- seed_start + iteration_id
  set.seed(current_seed)
}

# --- data gen -----------------------------------------------------------------
iter_start <- Sys.time()
cat("Generating data...\n"); data_start <- Sys.time()

if (sim_type == "acic") {
  n <- 1000; p <- 10
  df    <- gen_df_acic(
    model.trt="step", root.trt=0.35, overlap.trt="full",
    model.rsp="step", alignment=0.75, te.hetero="high",
    random.seed=current_seed, n=n, p=p
  )
  covs  <- df %>% dplyr::select(starts_with("X")) %>% colnames()
  form  <- as.formula(paste0("Z ~ ", paste0(covs, collapse = "+")))
  nbins <- 5

} else if (sim_type == "hainmueller") {
  nc <- 250; nt <- 50
  df    <- gen_df_hain(nc=nc, nt=nt, sigma_e="n100", outcome="nl2", sigma_y=1, ATE=0)
  form  <- as.formula("Z ~ X1 + X2 + X3 + X4 + X5 + X6")
  nbins <- 5

} else if (sim_type == "kang") {
  n <- 1000
  df    <- gen_df_kang(n=n)
  covs  <- df %>% dplyr::select(starts_with("X")) %>% colnames()
  form  <- as.formula(paste0("Z ~ ", paste0(covs, collapse = "+")))
  nbins <- 5

} else if (sim_type == "toy") {
  nc <- 500; nt <- 100; f0_sd <- 0.5
  df <- gen_df_adv(
    nc = nc, nt = nt, f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) { 3*X1 + 3*X2 },
    f0_fun = function(x, y) {
      matrix(c(x,y), ncol=2) %>%
        dmvnorm(mean=c(0.5,0.5), sigma=matrix(c(1,0.8,0.8,1),2)) * 20
    }
  )
  form  <- as.formula("Z ~ X1 + X2")
  nbins <- 6

} else {
  stop("sim_type must be one of: 'acic', 'hainmueller', 'kang', 'toy'")
}

# --- scaling + feasibility filter ---------------------------------------------
covs <- parse_form(form)$covs
dist_scaling <- df %>%
  summarize(across(
    all_of(covs),
    function(x) if (is.numeric(x)) nbins / (max(x) - min(x)) else 1000
  ))
preds_csm <- get_cal_matches(data=df, metric="maximum", scaling=dist_scaling,
                             rad_method="fixed", est_method="average", k=25)
preds_cem <- get_cem_matches(data=df, num_bins=nbins, est_method="average", return="all")
ninf      <- length(attr(preds_csm, "unmatched_units"))
ninf_cem  <- sum(df$Z) - length(attr(preds_cem, "feasible_units"))
# df        <- df %>% filter(!id %in% attr(preds_csm, "unmatched_units"))
data_time <- as.numeric(difftime(Sys.time(), data_start, units = "secs"))

true_ATT <- if (sim_type %in% c("acic","toy")) {
  df %>% dplyr::filter(Z == 1) %>% summarize(att = mean(Y1 - Y0)) %>% pull(att)
} else 0

cat(sprintf("Data generated: n=%d, true_ATT=%.3f, time=%.2fs\n", nrow(df), true_ATT, data_time))

# --- run methods --------------------------------------------------------------
method_results <- run_all_methods(
  df               = df,
  form             = form,
  dist_scaling     = dist_scaling,
  nbins            = nbins,
  selected_methods = SELECTED_METHODS
)

# --- save ---------------------------------------------------------------------
total_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

res <- tibble(
  runid             = iteration_id,
  seed              = current_seed,
  sim_type          = sim_type,
  elapsed_time_secs = total_time,
  ninf              = ninf,
  ninf_cem          = ninf_cem,
  true_ATT          = true_ATT
) %>%
  bind_cols(
    method_results %>%
      dplyr::select(method, ATT_est) %>%
      tidyr::pivot_wider(names_from = method, values_from = ATT_est)
  )

output_file <- file.path(output_dir, sprintf("iter_%04d.csv", iteration_id))

if (file.exists(output_file)) {
  existing <- read_csv(output_file, show_col_types = FALSE)
  cat(sprintf("Merging into existing file (updating: %s)\n", paste(SELECTED_METHODS, collapse = ", ")))
  for (col in SELECTED_METHODS) existing[[col]] <- res[[col]]
  res <- existing
}

write_csv(res, output_file)
cat(sprintf("\nIteration %d complete! Total time: %.2f seconds\n", iteration_id, total_time))
cat(sprintf("Results saved to: %s\n", output_file))
