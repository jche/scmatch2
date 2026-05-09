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
  "or_lm", "or_lm_main", "or_bart",
  "ps_lm", "ps_bart",
  "csm_scm", "csm_avg", "cem_scm", "cem_avg", "onenn",
  "csm_scm_half", "csm_avg_half", "cem_scm_half", "cem_avg_half",
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
      matching_half  = c("csm_scm_half", "csm_avg_half", "cem_scm_half", "cem_avg_half"),
      ps             = c("ps_lm", "ps_bart"),
      or             = c("or_lm", "or_lm_main", "or_bart"),
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
#' @param form       Formula of the form \code{Z ~ X1 + X2 + ...} naming the
#'                   treatment (LHS) and all covariates (RHS). All sub-formulas
#'                   required by individual methods are derived from this:
#'                   \itemize{
#'                     \item PS main-effects (bal1, twang): \code{form} as-is.
#'                     \item PS with interactions (bal2, ps_lm): \code{update(form, . ~ (.)^2)}.
#'                     \item Outcome model (or_lm): \code{update(form, Y ~ (.)^2)}.
#'                     \item Covariate matrix (or_bart, ps_bart, tmle, aipw, cf, kbal):
#'                           RHS variables extracted via \code{parse_form(form)$covs}.
#'                   }
#' @param dist_scaling Named numeric vector/row giving per-covariate scaling
#'                   (output of the summarize(across(...)) call in
#'                   run_single_iteration.R).
#' @param nbins      Number of CEM bins (used for cem_scm / cem_avg).
#' @param selected_methods Character vector of method names or group shorthands
#'                   to run (see \code{expand_method_groups}). Defaults to all
#'                   fast methods.
#' @param verbose    Logical. If \code{TRUE} (default), print method names and
#'                   estimates as they are computed. Set to \code{FALSE} to
#'                   suppress all output (useful in tests or non-interactive
#'                   scripts).
#'
#' @return Tibble with columns: method, ATT_est, secs.
run_all_methods <- function(df,
                            form,
                            dist_scaling,
                            nbins,
                            selected_methods = "fast",
                            verbose = TRUE) {

  selected_methods <- expand_method_groups(selected_methods)

  # Derive all sub-formulas from the single treatment formula
  covs      <- parse_form(form)$covs
  form_int  <- update(form, . ~ (.)^2)   # PS formula with all pairwise interactions
  yform_int <- update(form, Y ~ (.)^2)   # outcome formula with interactions (for or_lm)
  yform_main <- update(form, Y ~ .)      # outcome formula with main effects only (for or_lm_main)

  safe_compute <- function(method_name, expr) {
    if (!(method_name %in% selected_methods)) {
      return(list(est = NA_real_, secs = NA_real_))
    }
    if (verbose)
      cat(sprintf("  %-20s ... ", method_name))
    t0 <- Sys.time()

    result <- tryCatch(expr, error = function(e) {
      if (verbose)
        cat(sprintf("ERROR: %s\n", e$message))
      NA_real_
    })

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (verbose && !is.na(result))
      cat(sprintf("%.4f  (%.1fs)\n", result, elapsed))
    list(est = result, secs = elapsed)
  }

  if (verbose) cat("Computing estimates:\n")

  r_diff    <- safe_compute("diff",    get_att_diff(df))
  r_bal1    <- safe_compute("bal1",    get_att_bal(df, form,     rep(0.01, length(covs))))
  r_bal2    <- safe_compute("bal2",    get_att_bal(df, form_int, rep(0.1,  length(covs) + choose(length(covs), 2))))
  r_or_lm      <- safe_compute("or_lm",      get_att_or_lm(df, form = yform_int))
  r_or_lm_main <- safe_compute("or_lm_main", get_att_or_lm(df, form = yform_main))
  r_or_bart <- safe_compute("or_bart", get_att_or_bart(df, form = form))
  r_ps_lm   <- safe_compute("ps_lm",  get_att_ps_lm(df,  form_int))
  r_ps_bart <- safe_compute("ps_bart", get_att_ps_bart(df, form = form))
  # Fair comparison (Proposition §sec:compCEM): CSM caliper = range/(2*nbins),
  # CEM bin-width = range/nbins = 2 × CSM caliper → equal catchment, CSM bias 2× tighter.
  dist_scaling_csm <- df %>%
    dplyr::summarize(dplyr::across(
      dplyr::all_of(covs),
      function(x) if (is.numeric(x)) (2 * nbins) / (max(x) - min(x)) else 1000
    ))
  r_csm_scm <- safe_compute("csm_scm", get_att_csm(df, scaling = dist_scaling_csm, est_method = "scm",     rad_method = "fixed"))
  r_csm_avg <- safe_compute("csm_avg", get_att_csm(df, scaling = dist_scaling_csm, est_method = "average", rad_method = "fixed"))
  r_cem_scm <- safe_compute("cem_scm", get_att_cem(df, num_bins = nbins, est_method = "scm",     estimand = "CEM-ATT"))
  r_cem_avg <- safe_compute("cem_avg", get_att_cem(df, num_bins = nbins, est_method = "average", estimand = "CEM-ATT"))
  r_onenn   <- safe_compute("onenn",   get_att_1nn(df, scaling = dist_scaling))

  # Half-bin variants at half granularity: CEM uses nbins/2, CSM uses 2*(nbins/2) = nbins
  # (dist_scaling is already nbins/(max-min), so it gives CSM caliper = range/nbins = 2× CEM bin-width)
  r_csm_scm_half <- safe_compute("csm_scm_half", get_att_csm(df, scaling = dist_scaling,      est_method = "scm",     rad_method = "fixed"))
  r_csm_avg_half <- safe_compute("csm_avg_half", get_att_csm(df, scaling = dist_scaling,      est_method = "average", rad_method = "fixed"))
  r_cem_scm_half <- safe_compute("cem_scm_half", get_att_cem(df, num_bins = nbins / 2, est_method = "scm",     estimand = "CEM-ATT"))
  r_cem_avg_half <- safe_compute("cem_avg_half", get_att_cem(df, num_bins = nbins / 2, est_method = "average", estimand = "CEM-ATT"))
  r_tmle1   <- safe_compute("tmle1",   get_att_tmle(df, form = form, Q.SL.library = SL.library1,  g.SL.library = SL.library1))
  r_aipw1   <- safe_compute("aipw1",   get_att_aipw(df, form = form, Q.SL.library = SL.library1,  g.SL.library = SL.library1))
  r_tmle2   <- safe_compute("tmle2",   get_att_tmle(df, form = form, Q.SL.library = SL.library3Q, g.SL.library = SL.library3g))
  r_aipw2   <- safe_compute("aipw2",   get_att_aipw(df, form = form, Q.SL.library = SL.library2,  g.SL.library = SL.library2))
  r_cf      <- safe_compute("causal_forest", get_att_causal_forest(df, form = form))
  r_twang   <- safe_compute("twang",   get_att_twang(df, form = form))
  r_kbal    <- safe_compute("kbal",    get_att_kbal(df,  form = form))

  tibble::tibble(
    method  = ALL_METHODS,
    ATT_est = c(r_diff$est, r_bal1$est, r_bal2$est,
                r_or_lm$est, r_or_lm_main$est, r_or_bart$est,
                r_ps_lm$est, r_ps_bart$est,
                r_csm_scm$est, r_csm_avg$est, r_cem_scm$est, r_cem_avg$est, r_onenn$est,
                r_csm_scm_half$est, r_csm_avg_half$est, r_cem_scm_half$est, r_cem_avg_half$est,
                r_tmle1$est, r_aipw1$est, r_tmle2$est, r_aipw2$est,
                r_cf$est, r_twang$est, r_kbal$est),
    secs    = c(r_diff$secs, r_bal1$secs, r_bal2$secs,
                r_or_lm$secs, r_or_lm_main$secs, r_or_bart$secs,
                r_ps_lm$secs, r_ps_bart$secs,
                r_csm_scm$secs, r_csm_avg$secs, r_cem_scm$secs, r_cem_avg$secs, r_onenn$secs,
                r_csm_scm_half$secs, r_csm_avg_half$secs, r_cem_scm_half$secs, r_cem_avg_half$secs,
                r_tmle1$secs, r_aipw1$secs, r_tmle2$secs, r_aipw2$secs,
                r_cf$secs, r_twang$secs, r_kbal$secs)
  )
}
