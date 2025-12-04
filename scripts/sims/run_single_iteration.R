#!/usr/bin/env Rscript
# Run a single iteration of the comprehensive simulation (+ toy)
# Usage:
# Rscript scripts/sims/canonical_run_single_iteration.R <sim_type> <iteration_id> [methods...]
# sim_type: one of 'acic','hainmueller','kang','toy'
# iteration_id: integer
# methods: optional list and/or comma-separated names (see supported list below)

devtools::load_all()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript canonical_run_single_iteration.R <sim_type> <iteration_id> [methods...]")
}
sim_type    <- args[1]
iteration_id <- as.integer(args[2])

# --- parse optional methods ---------------------------------------------------
raw_methods <- character(0)
if (length(args) > 2) {
  margs <- args[-(1:2)]
  # allow comma-separated and space-separated
  raw_methods <- unlist(strsplit(margs, split = ","))
  raw_methods <- trimws(tolower(raw_methods))
  raw_methods <- raw_methods[nzchar(raw_methods)]
}

# canonical method list
ALL_METHODS <- c(
  "diff","bal1","bal2",
  "or_lm","or_bart",
  "ps_lm","ps_bart",
  "csm_scm","csm_avg","cem_scm","cem_avg","onenn",
  "tmle1","aipw1","tmle2","aipw2",
  "causal_forest","twang","kbal"
)

# expandable groups
expand_group <- function(keys) {
  out <- c()
  for (k in keys) {
    if (k %in% c("all")) {
      out <- c(out, ALL_METHODS)
    } else if (k %in% c("old")) {
      out <- c(out, setdiff(ALL_METHODS, c("causal_forest","twang","kbal")))
    } else if (k %in% c("new")) {
      out <- c(out, c("causal_forest","twang","kbal"))
    } else if (k %in% c("matching")) {
      out <- c(out, c("csm_scm","csm_avg","cem_scm","cem_avg","onenn"))
    } else if (k %in% c("ps")) {
      out <- c(out, c("ps_lm","ps_bart"))
    } else if (k %in% c("or")) {
      out <- c(out, c("or_lm","or_bart"))
    } else if (k %in% c("bal","balance")) {
      out <- c(out, c("bal1","bal2"))
    } else if (k %in% c("dr","doubly_robust")) {
      out <- c(out, c("tmle1","aipw1","tmle2","aipw2"))
    } else if (k %in% ALL_METHODS) {
      out <- c(out, k)
    } else {
      warning(sprintf("Unknown method/group '%s' ignored.", k))
    }
  }
  unique(out)
}

SELECTED_METHODS <- if (length(raw_methods) == 0) ALL_METHODS else expand_group(raw_methods)

cat(sprintf("\n=== Running %s simulation, iteration %d ===\n", sim_type, iteration_id))
if (length(raw_methods) == 0) {
  cat("Methods: ALL\n")
} else {
  cat("Methods: ", paste(SELECTED_METHODS, collapse = ", "), "\n")
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
source(here::here("scripts/wrappers.R"))
source(here::here("R/utils.R"))
source(here::here("R/sim_data.R"))

SL.library1  <- c("SL.mean", "SL.lm", "SL.glm")
SL.library2  <- c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.xgboost")
SL.library3Q <- c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")
SL.library3g <- c("SL.glm", "tmle.SL.dbarts.k.5", "SL.gam")

output_dir <- file.path("data", "outputs", "sims", sim_type)
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
  df <- gen_df_acic(
    model.trt="step", root.trt=0.35, overlap.trt="full",
    model.rsp="step", alignment=0.75, te.hetero="high",
    random.seed=current_seed, n=n, p=p
  )
  covs    <- df %>% select(starts_with("X")) %>% colnames()
  zform1  <- as.formula(paste0("Z ~ ", paste0(covs, collapse = "+")))
  zform2  <- as.formula(paste0("Z ~ (", paste0(covs, collapse = "+"), ")^2"))
  form2   <- as.formula(paste0("Y ~ (", paste0(covs, collapse = "+"), ")^2"))
  nbins   <- 5

} else if (sim_type == "hainmueller") {
  nc <- 250; nt <- 50
  df <- gen_df_hain(nc=nc, nt=nt, sigma_e="n100", outcome="nl2", sigma_y=1, ATE=0)
  covs   <- c("X1","X2","X3","X4","X5","X6")
  zform1 <- as.formula("Z ~ X1+X2+X3+X4+X5+X6")
  zform2 <- as.formula("Z ~ (X1+X2+X3+X4+X5+X6)^2")
  form2  <- as.formula("Y ~ (X1+X2+X3+X4+X5+X6)^2")
  nbins  <- 5

} else if (sim_type == "kang") {
  n <- 1000
  df <- gen_df_kang(n=n)
  covs   <- df %>% select(starts_with("X")) %>% colnames()
  zform1 <- as.formula(paste0("Z ~ ", paste0(covs, collapse = "+")))
  zform2 <- as.formula(paste0("Z ~ (", paste0(covs, collapse = "+"), ")^2"))
  form2  <- as.formula(paste0("Y ~ (", paste0(covs, collapse = "+"), ")^2"))
  nbins  <- 5

} else if (sim_type == "toy") {
  nc <- 500; nt <- 100; f0_sd <- 0.5
  df <- gen_df_adv(
    nc = nc, nt = nt, f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) {3*X1 + 3*X2},
    f0_fun = function(x, y) {
      matrix(c(x,y), ncol=2) %>%
        dmvnorm(mean=c(0.5,0.5), sigma=matrix(c(1,0.8,0.8,1),2)) * 20
    }
  )
  covs   <- c("X1","X2")
  zform1 <- as.formula("Z ~ X1 + X2")
  zform2 <- as.formula("Z ~ X1 + X2")
  form2  <- as.formula("Y ~ X1 + X2")
  nbins  <- 6
} else {
  stop("sim_type must be one of: 'acic', 'hainmueller', 'kang', 'toy'")
}

# common scaling + CSM/CEM filtering
dist_scaling <- df %>%
  summarize(across(
    starts_with("X"),
    function(x) if (is.numeric(x)) (if (sim_type=="toy") 2*nbins else nbins) / (max(x)-min(x)) else 1000
  ))
preds_csm <- get_cal_matches(data=df, metric="maximum", scaling=dist_scaling,
                             rad_method="fixed", est_method="average", k=25)
preds_cem <- get_cem_matches(data=df, num_bins=nbins, est_method="average", return="all")
ninf      <- length(attr(preds_csm, "unmatched_units"))
ninf_cem  <- sum(df$Z) - length(attr(preds_cem, "feasible_units"))
df <- df %>% filter(!id %in% attr(preds_csm, "unmatched_units"))
data_time <- as.numeric(difftime(Sys.time(), data_start, units = "secs"))

true_ATT <- if (sim_type %in% c("acic","toy")) {
  df %>% dplyr::filter(Z) %>% summarize(att = mean(Y1 - Y0)) %>% pull(att)
} else 0

cat(sprintf("Data generated: n=%d, true_ATT=%.3f, time=%.2fs\n", nrow(df), true_ATT, data_time))

# --- compute only selected methods -------------------------------------------
safe_compute <- function(method_name, expr) {
  t_start <- Sys.time()
  result <- tryCatch(expr, error = function(e) {
    cat(sprintf("  %s: ERROR - %s\n", method_name, e$message)); return(NA_real_)
  })
  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (!is.na(result)) cat(sprintf("  %s: %.3f (%.2fs)\n", method_name, result, t_elapsed))
  result
}
run_if <- function(name, code) if (name %in% SELECTED_METHODS) safe_compute(name, code) else NA_real_

cat("\nComputing estimates...\n")
att_diff     <- run_if("diff",    get_att_diff(df))
att_bal1     <- run_if("bal1",    get_att_bal(df, zform1, rep(0.01, length(covs))))
att_bal2     <- run_if("bal2",    get_att_bal(df, zform2, rep(0.1,  length(covs) + choose(length(covs),2))))
att_or_lm    <- run_if("or_lm",   get_att_or_lm(df, form=form2))
att_or_bart  <- run_if("or_bart", get_att_or_bart(df, covs=covs))
att_ps_lm    <- run_if("ps_lm",   get_att_ps_lm(df, zform2))
att_ps_bart  <- run_if("ps_bart", get_att_ps_bart(df, covs=covs))
att_csm_scm  <- run_if("csm_scm", get_att_csm(df, scaling=dist_scaling, est_method="scm",     rad_method="fixed"))
att_csm_avg  <- run_if("csm_avg", get_att_csm(df, scaling=dist_scaling, est_method="average", rad_method="fixed"))
att_cem_scm  <- run_if("cem_scm", get_att_cem(df, num_bins=nbins, est_method="scm", estimand="CEM-ATT"))
att_cem_avg  <- run_if("cem_avg", get_att_cem(df, num_bins=nbins, est_method="average", estimand="CEM-ATT"))
att_onenn    <- run_if("onenn",   get_att_1nn(df, scaling=dist_scaling))
att_tmle1    <- run_if("tmle1",   get_att_tmle(df, covs=covs, Q.SL.library=SL.library1,  g.SL.library=SL.library1))
att_aipw1    <- run_if("aipw1",   get_att_aipw(df, covs=covs, Q.SL.library=SL.library1,  g.SL.library=SL.library1))
att_tmle2    <- run_if("tmle2",   get_att_tmle(df, covs=covs, Q.SL.library=SL.library3Q, g.SL.library=SL.library3g))
att_aipw2    <- run_if("aipw2",   get_att_aipw(df, covs=covs, Q.SL.library=SL.library2,  g.SL.library=SL.library2))
att_cf       <- run_if("causal_forest", get_att_causal_forest(df, covs=covs))
att_twang    <- run_if("twang",        get_att_twang(df, form=zform1))
att_kbal     <- run_if("kbal",         get_att_kbal(df,  covs=covs))

# --- save ---------------------------------------------------------------------
total_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
res <- tibble(
  runid = iteration_id, seed = current_seed, sim_type = sim_type,
  elapsed_time_secs = total_time, ninf = ninf, ninf_cem = ninf_cem, true_ATT = true_ATT,
  diff = att_diff, bal1 = att_bal1, bal2 = att_bal2, or_lm = att_or_lm, or_bart = att_or_bart,
  ps_lm = att_ps_lm, ps_bart = att_ps_bart, csm_scm = att_csm_scm, csm_avg = att_csm_avg,
  cem_scm = att_cem_scm, cem_avg = att_cem_avg, onenn = att_onenn, tmle1 = att_tmle1,
  aipw1 = att_aipw1, tmle2 = att_tmle2, aipw2 = att_aipw2, causal_forest = att_cf,
  twang = att_twang, kbal = att_kbal
)
output_file <- file.path(output_dir, sprintf("iter_%04d.csv", iteration_id))
write_csv(res, output_file)
cat(sprintf("\nIteration %d complete! Total time: %.2f seconds\n", iteration_id, total_time))
cat(sprintf("Results saved to: %s\n", output_file))
