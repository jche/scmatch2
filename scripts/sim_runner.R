# scripts/sim_runner.R
#
# Orchestration layer for the bias/MSE simulation study.
# Provides:
#   - ALL_METHODS        : canonical ordered character vector of method names
#   - expand_method_groups() : expand shorthand group names to method lists
#   - run_all_methods()      : run a set of methods on a prepped dataset,
#                              returning a tidy tibble of estimates + timings
#
# Low-level get_att_*() wrappers live in scripts/wrappers.R.
# Data-generation functions live in R/sim_data.R (loaded via devtools::load_all()).

# ------------------------------------------------------------------------------
# Method registry
# ------------------------------------------------------------------------------

ALL_METHODS <- c(
  "diff", "bal1", "bal2",
  "or_lm", "or_bart",
  "ps_lm", "ps_bart",
  "csm_scm", "csm_avg", "cem_scm", "cem_avg", "onenn",
  "tmle1", "aipw1", "tmle2", "aipw2",
  "causal_forest", "twang", "kbal"
)

SLOW_METHODS <- c("tmle1", "aipw1", "tmle2", "aipw2", "causal_forest", "twang", "kbal")

#' Expand shorthand group names into a concrete vector of method names.
#'
#' @param keys Character vector of method names and/or group shorthands.
#'   Recognized groups: "all", "old", "new", "matching", "ps", "or", "bal",
#'   "balance", "dr", "doubly_robust", "slow", "fast".
#' @return Character vector of unique, validated method names.
expand_method_groups <- function(keys, drop_slow = FALSE) {
  out <- character(0)
  for (k in tolower(trimws(keys))) {
    out <- c(out, switch(k,
      all            = ALL_METHODS,
      old            = setdiff(ALL_METHODS, c("causal_forest", "twang", "kbal")),
      new            = c("causal_forest", "twang", "kbal"),
      matching       = c("csm_scm", "csm_avg", "cem_scm", "cem_avg", "onenn"),
      ps             = c("ps_lm", "ps_bart"),
      or             = c("or_lm", "or_bart"),
      bal            = ,
      balance        = c("bal1", "bal2"),
      dr             = ,
      doubly_robust  = c("tmle1", "aipw1", "tmle2", "aipw2"),
      slow           = SLOW_METHODS,
      fast           = setdiff(ALL_METHODS, SLOW_METHODS),
      {
        if (k %in% ALL_METHODS) k
        else { warning(sprintf("Unknown method/group '%s' ignored.", k)); character(0) }
      }
    ))
  }

  if ( drop_slow ) {
    out <- setdiff(out, SLOW_METHODS)
  }

  unique(out)
}

# ------------------------------------------------------------------------------
# SuperLearner libraries (used by tmle / aipw methods)
# ------------------------------------------------------------------------------

SL.library1  <- c("SL.mean", "SL.lm", "SL.glm")
SL.library2  <- c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.xgboost")
SL.library3Q <- c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")
SL.library3g <- c("SL.glm", "tmle.SL.dbarts.k.5", "SL.gam")

# ------------------------------------------------------------------------------
# Core runner
# ------------------------------------------------------------------------------

#' Run ATT estimation methods on a prepared dataset.
#'
#' The dataset should already have infeasible treated units removed (i.e., units
#' with no CSM match). Call get_cal_matches() upstream and filter before passing
#' df here. See run_single_iteration.R for the standard preparation steps.
#'
#' @param df         Data frame with columns Y, Z, id, and covariates.
#' @param covs       Character vector of covariate column names (e.g. c("X1","X2")).
#' @param zform1     PS formula with main effects (e.g. Z ~ X1 + X2).
#' @param zform2     PS formula with interactions (e.g. Z ~ X1*X2).
#' @param form2      Outcome formula for or_lm (e.g. Y ~ X1 + X2).
#' @param dist_scaling Named numeric vector/row giving per-covariate scaling
#'                   (output of the summarize(across(...)) call in run_single_iteration.R).
#' @param nbins      Number of CEM bins (used for cem_scm / cem_avg).
#' @param selected_methods Character vector of method names to run. Defaults to
#'                   all fast methods. Pass expand_method_groups("all") for every
#'                   method including slow ones.
#'
#' @return Tibble with columns: method, ATT_est, secs.
run_all_methods <- function(df,
                            covs,
                            zform1,
                            zform2,
                            form2,
                            dist_scaling,
                            nbins,
                            selected_methods = "fast" ) {

  selected_methods <- expand_method_groups(selected_methods)

  safe_compute <- function(method_name, expr) {
    if (!(method_name %in% selected_methods)) {
      return(list(est = NA_real_, secs = NA_real_))
    }
    cat(sprintf("  %-20s ... ", method_name))
    t0 <- Sys.time()
    result <- tryCatch(expr, error = function(e) {
      cat(sprintf("ERROR: %s\n", e$message))
      NA_real_
    })
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.na(result)) cat(sprintf("%.4f  (%.1fs)\n", result, elapsed))
    list(est = result, secs = elapsed)
  }

  cat("Computing estimates:\n")

  r_diff    <- safe_compute("diff",    get_att_diff(df))
  r_bal1    <- safe_compute("bal1",    get_att_bal(df, zform1, rep(0.01, length(covs))))
  r_bal2    <- safe_compute("bal2",    get_att_bal(df, zform2, rep(0.1, length(covs) + choose(length(covs), 2))))
  r_or_lm   <- safe_compute("or_lm",  get_att_or_lm(df, form = form2))
  r_or_bart <- safe_compute("or_bart", get_att_or_bart(df, covs = covs))
  r_ps_lm   <- safe_compute("ps_lm",  get_att_ps_lm(df, zform2))
  r_ps_bart <- safe_compute("ps_bart", get_att_ps_bart(df, covs = covs))
  r_csm_scm <- safe_compute("csm_scm", get_att_csm(df, scaling = dist_scaling, est_method = "scm",     rad_method = "fixed"))
  r_csm_avg <- safe_compute("csm_avg", get_att_csm(df, scaling = dist_scaling, est_method = "average", rad_method = "fixed"))
  r_cem_scm <- safe_compute("cem_scm", get_att_cem(df, num_bins = nbins, est_method = "scm",     estimand = "CEM-ATT"))
  r_cem_avg <- safe_compute("cem_avg", get_att_cem(df, num_bins = nbins, est_method = "average", estimand = "CEM-ATT"))
  r_onenn   <- safe_compute("onenn",   get_att_1nn(df, scaling = dist_scaling))
  r_tmle1   <- safe_compute("tmle1",   get_att_tmle(df, covs = covs, Q.SL.library = SL.library1,  g.SL.library = SL.library1))
  r_aipw1   <- safe_compute("aipw1",   get_att_aipw(df, covs = covs, Q.SL.library = SL.library1,  g.SL.library = SL.library1))
  r_tmle2   <- safe_compute("tmle2",   get_att_tmle(df, covs = covs, Q.SL.library = SL.library3Q, g.SL.library = SL.library3g))
  r_aipw2   <- safe_compute("aipw2",   get_att_aipw(df, covs = covs, Q.SL.library = SL.library2,  g.SL.library = SL.library2))
  r_cf      <- safe_compute("causal_forest", get_att_causal_forest(df, covs = covs))
  r_twang   <- safe_compute("twang",   get_att_twang(df, form = zform1))
  r_kbal    <- safe_compute("kbal",    get_att_kbal(df, covs = covs))

  tibble::tibble(
    method  = ALL_METHODS,
    ATT_est = c(r_diff$est, r_bal1$est, r_bal2$est,
                r_or_lm$est, r_or_bart$est,
                r_ps_lm$est, r_ps_bart$est,
                r_csm_scm$est, r_csm_avg$est, r_cem_scm$est, r_cem_avg$est, r_onenn$est,
                r_tmle1$est, r_aipw1$est, r_tmle2$est, r_aipw2$est,
                r_cf$est, r_twang$est, r_kbal$est),
    secs    = c(r_diff$secs, r_bal1$secs, r_bal2$secs,
                r_or_lm$secs, r_or_bart$secs,
                r_ps_lm$secs, r_ps_bart$secs,
                r_csm_scm$secs, r_csm_avg$secs, r_cem_scm$secs, r_cem_avg$secs, r_onenn$secs,
                r_tmle1$secs, r_aipw1$secs, r_tmle2$secs, r_aipw2$secs,
                r_cf$secs, r_twang$secs, r_kbal$secs)
  )
}
