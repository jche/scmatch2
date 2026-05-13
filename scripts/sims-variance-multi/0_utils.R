# scripts/sims-variance-multi/0_utils.R
#
# Utilities for the 3-factor variance-estimator simulation study.
#
# Factors
# -------
#   1. Overlap (5 levels):        very_low, low, mid, high, very_high
#                                 controlled via prop_nc_unif
#   2. Error structure (2 levels): homo  — σ₀(x) = 0.5  (constant)
#                                  het   — σ₀(x) bimodal Gaussian bumps
#   3. Common variance (2 levels): common    — σ₁ = σ₀  (sigma1_extra = 0)
#                                  no_common — σ₁ = σ₀ + 2 (sigma1_extra = 2)
#
# 5 × 2 × 2 = 20 design cells, 998 replications.
#
# Four variance estimators
# ------------------------
#   homo       get_ATT_estimate(variance_method = "pooled")
#              SE = sqrt(V_E + V_P)   [total variance]
#
#   het        get_ATT_estimate(variance_method = "pooled_het")
#              SE = sqrt(V_E + V_P)   [total variance]
#
#   alt_common get_measurement_error_variance_alt(use_common_variance = TRUE)
#              SE = sqrt(V_E_alt)     [measurement error only]
#
#   alt_tt     get_measurement_error_variance_alt(use_common_variance = FALSE, K = 2)
#              SE = sqrt(V_E_alt)     [measurement error only]
#
# Performance measure: coverage of τ_SATT
#   covered = (CI_lower ≤ SATT) & (SATT ≤ CI_upper)

devtools::load_all()

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
  library(here)
  library(mvtnorm)
})

# ─────────────────────────────────────────────────────────────────────────────
# Design constants
# ─────────────────────────────────────────────────────────────────────────────

OVERLAP_LABELS <- c("very_low", "low", "mid", "high", "very_high")

PROP_NC_UNIF <- c(
  very_low  = 1/10,
  low       = 1/5,
  mid       = 1/3,
  high      = 1/2,
  very_high = 2/3
)

DESIGN_GRID <- expand_grid(
  overlap_label = OVERLAP_LABELS,
  error_type    = c("homo", "het"),
  sigma1_extra  = c(0, 2)
) %>%
  mutate(
    common_label = if_else(sigma1_extra == 0, "common", "no_common"),
    # factor ordering for plots
    overlap_label = factor(overlap_label, levels = OVERLAP_LABELS),
    error_type    = factor(error_type,    levels = c("homo", "het")),
    common_label  = factor(common_label,  levels = c("common", "no_common"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# σ₀(x): bimodal Gaussian bumps (v2 parameters, A=5, h_sq=0.08, sigma_min=0.2)
# At treated cluster centres (0.25,0.25) and (0.75,0.75): σ₀ ≈ 5.2
# At control cluster centres: σ₀ ≈ 0.2
# ─────────────────────────────────────────────────────────────────────────────

sigma0_bimodal <- function(x1, x2,
                            sigma_min = 0.2, A = 5.0, h_sq = 0.08) {
  d1 <- (x1 - 0.25)^2 + (x2 - 0.25)^2
  d2 <- (x1 - 0.75)^2 + (x2 - 0.75)^2
  sigma_min + A * (exp(-d1 / h_sq) + exp(-d2 / h_sq))
}

# ─────────────────────────────────────────────────────────────────────────────
# Data generation
# ─────────────────────────────────────────────────────────────────────────────

#' Generate one dataset for the multifactor simulation.
#'
#' Uses gen_df_adv_k() with f0_sd = 0 to get the noiseless skeleton
#' (Y0_denoised, Y1_denoised), then adds noise according to error_type
#' and sigma1_extra.
#'
#' @param nc,nt  Control and treated sample sizes.
#' @param prop_nc_unif  Fraction of control units drawn uniformly
#'   (controls degree of overlap).
#' @param ctr_dist  Cluster separation parameter.
#' @param error_type  "homo" (σ₀ = 0.5 constant) or "het" (bimodal σ₀).
#' @param sigma1_extra  Extra treatment-side SD: σ₁ = σ₀ + sigma1_extra.
#'   0 = common variance; 2 = no-common variance.
#' @param seed  Optional seed (set before gen_df_adv_k call).
#' @return Data frame with Y, Y0, Y1, Z, X1, X2, id, sigma0, sigma1.
make_df_multi <- function(
    nc = 500, nt = 100,
    prop_nc_unif,
    ctr_dist     = 0.5,
    error_type   = c("homo", "het"),
    sigma1_extra = 0,
    seed         = NULL) {

  error_type <- match.arg(error_type)
  if (!is.null(seed)) set.seed(seed)

  f0_fun_mat <- function(X) {
    X <- as.matrix(X)
    mvtnorm::dmvnorm(X, mean = c(0.5, 0.5),
                     sigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)) * 20
  }
  tx_fun_mat <- function(X) {
    X <- as.matrix(X)
    3 * X[, 1] + 3 * X[, 2]
  }

  # Generate noiseless skeleton (f0_sd = 0 → Y0_denoised = f0_fun(X))
  df_raw <- gen_df_adv_k(
    nc = nc, nt = nt, k = 2,
    f0_sd         = 0,
    tx_effect_fun = tx_fun_mat,
    f0_fun        = f0_fun_mat,
    ctr_dist      = ctr_dist,
    prop_nc_unif  = prop_nc_unif
  )

  # σ₀ shape depends on error_type factor
  sig0 <- if (error_type == "homo") {
    rep(0.5, nrow(df_raw))
  } else {
    sigma0_bimodal(df_raw$X1, df_raw$X2)
  }

  sig1 <- sig0 + sigma1_extra

  eps0 <- rnorm(nrow(df_raw), 0, sig0)
  # When common variance (sigma1_extra = 0): reuse eps0 so σ₁ = σ₀ exactly.
  # When no-common (sigma1_extra > 0): draw independent treatment noise.
  eps1 <- if (sigma1_extra == 0) eps0 else rnorm(nrow(df_raw), 0, sig1)

  df_raw %>%
    mutate(
      Y0     = Y0_denoised + eps0,
      Y1     = Y1_denoised + eps1,
      Y      = ifelse(Z, Y1, Y0),
      Z      = as.integer(Z),
      sigma0 = sig0,
      sigma1 = sig1
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# Distance scaling  (2 * nbins / range, matching sims-variance-het convention)
# ─────────────────────────────────────────────────────────────────────────────

compute_scaling <- function(df, nbins = 6) {
  df %>%
    summarise(across(
      starts_with("X"),
      function(x) if (is.numeric(x)) (2 * nbins) / (max(x) - min(x)) else 1000
    ))
}

# ─────────────────────────────────────────────────────────────────────────────
# One iteration: match + four variance estimators
# ─────────────────────────────────────────────────────────────────────────────

#' Run one simulation iteration for a single design cell.
#'
#' @param i         Replication index (for runID column).
#' @param overlap_label,error_type,sigma1_extra  Design cell identifiers.
#' @param prop_nc_unif  Looked up from PROP_NC_UNIF inside sim_master_multi;
#'   passed directly here.
#' @param seed_addition  Seed (unique per cell × replication).
#' @param K_tt  Treated-to-treated neighbours for alt_tt (default 2).
#' @param verbose  Print error messages from individual estimators.
#' @return Tibble with 4 rows (one per estimator) and all metadata columns.
one_iter <- function(
    i,
    overlap_label,
    error_type    = c("homo", "het"),
    sigma1_extra  = 0,
    prop_nc_unif,
    nc            = 500,
    nt            = 100,
    ctr_dist      = 0.5,
    nbins         = 6,
    k_match       = 2,
    K_tt          = 2,
    seed_addition,
    verbose       = FALSE
) {
  error_type <- match.arg(error_type)
  set.seed(seed_addition)

  df <- make_df_multi(
    nc = nc, nt = nt,
    prop_nc_unif = prop_nc_unif,
    ctr_dist     = ctr_dist,
    error_type   = error_type,
    sigma1_extra = sigma1_extra,
    seed         = seed_addition
  )

  scaling <- compute_scaling(df, nbins = nbins)

  mtch <- tryCatch(
    get_cal_matches(
      data       = df,
      Z ~ X1 + X2,
      rad_method = "adaptive",
      scaling    = scaling,
      k          = k_match,
      warn       = FALSE,
      est_method = "scm"
    ),
    error = function(e) {
      if (verbose) message("Matching failed iter ", i, ": ", e$message)
      NULL
    }
  )

  att_est   <- if (!is.null(mtch)) {
    tryCatch(get_att_point_est(mtch), error = function(e) NA_real_)
  } else NA_real_

  true_satt <- df %>%
    filter(Z == 1) %>%
    summarise(att = mean(Y1 - Y0)) %>%
    pull(att)

  # ── row constructors ────────────────────────────────────────────────────────
  make_row <- function(method_name, SE, V_E, V_P = NA_real_) {
    tibble(
      runID            = i,
      inference_method = method_name,
      att_est          = att_est,
      SE               = SE,
      CI_lower         = att_est - 1.96 * SE,
      CI_upper         = att_est + 1.96 * SE,
      V_E              = V_E,
      V_P              = V_P,
      SATT             = true_satt
    )
  }

  make_na_row <- function(method_name) {
    make_row(method_name, SE = NA_real_, V_E = NA_real_, V_P = NA_real_)
  }

  methods_all <- c("homo", "het", "alt_common", "alt_tt")

  if (is.null(mtch)) {
    return(
      bind_rows(lapply(methods_all, make_na_row)) %>%
        add_metadata(i, overlap_label, error_type, sigma1_extra, nc, nt, k_match, "Matching Failed")
    )
  }

  matches_full <- tryCatch(full_unit_table(mtch), error = function(e) NULL)

  # ── homo: homoskedastic pooled, SE = sqrt(V_E + V_P) ──────────────────────
  row_homo <- tryCatch({
    res <- get_ATT_estimate(mtch, variance_method = "pooled")
    make_row("homo", res$SE, res$V_E, res$V_P)
  }, error = function(e) {
    if (verbose) message("homo error: ", e$message)
    make_na_row("homo")
  })

  # ── het: heteroskedastic pooled, SE = sqrt(V_E + V_P) ─────────────────────
  row_het <- tryCatch({
    res <- get_ATT_estimate(mtch, variance_method = "pooled_het")
    make_row("het", res$SE, res$V_E, res$V_P)
  }, error = function(e) {
    if (verbose) message("het error: ", e$message)
    make_na_row("het")
  })

  # ── alt_common: control-side cluster variances, SE = sqrt(V_E_alt) ─────────
  row_alt_common <- tryCatch({
    if (is.null(matches_full)) stop("full_unit_table unavailable")
    alt <- get_measurement_error_variance_alt(matches_full,
                                              use_common_variance = TRUE)
    make_row("alt_common", sqrt(pmax(0, alt$V_E_alt)), alt$V_E_alt)
  }, error = function(e) {
    if (verbose) message("alt_common error: ", e$message)
    make_na_row("alt_common")
  })

  # ── alt_tt: treated-to-treated K-NN, SE = sqrt(V_E_alt) ───────────────────
  row_alt_tt <- tryCatch({
    if (is.null(matches_full)) stop("full_unit_table unavailable")
    alt_tt <- get_measurement_error_variance_alt(
      matches_table       = matches_full,
      df                  = df,
      use_common_variance = FALSE,
      K                   = K_tt,
      covs                = c("X1", "X2"),
      scaling             = scaling
    )
    make_row("alt_tt", sqrt(pmax(0, alt_tt$V_E_alt)), alt_tt$V_E_alt)
  }, error = function(e) {
    if (verbose) message("alt_tt error: ", e$message)
    make_na_row("alt_tt")
  })

  sz <- tryCatch(
    {
      s <- calc_N_T_N_C(full_unit_table(mtch))
      list(N_T = s$N_T, N_C_tilde = s$N_C_tilde)
    },
    error = function(e) list(N_T = NA_real_, N_C_tilde = NA_real_)
  )

  bind_rows(row_homo, row_het, row_alt_common, row_alt_tt) %>%
    mutate(N_T = sz$N_T, ESS_C = sz$N_C_tilde) %>%
    add_metadata(i, overlap_label, error_type, sigma1_extra, nc, nt, k_match, "Success")
}

# ── helper to attach design-cell metadata ───────────────────────────────────
add_metadata <- function(df, i, overlap_label, error_type, sigma1_extra,
                         nc, nt, k_match, status) {
  df %>%
    mutate(
      runID        = i,
      overlap_label = overlap_label,
      error_type   = error_type,
      sigma1_extra = sigma1_extra,
      common_label = if_else(sigma1_extra == 0, "common", "no_common"),
      nc           = nc,
      nt           = nt,
      k_match      = k_match,
      status       = status
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# sim_master_multi: one SLURM job = one iteration × all 20 design cells
# ─────────────────────────────────────────────────────────────────────────────

#' Run replication `iteration` across every row of DESIGN_GRID.
#'
#' Seeds are constructed to be unique per (iteration × cell) while remaining
#' reproducible.  The formula is:
#'   seed = iteration
#'        + 1000  * (overlap_idx - 1)
#'        + 10000 * (error_idx   - 1)   [1=homo, 2=het]
#'        + 50000 * (sigma1_idx  - 1)   [1=common, 2=no_common]
#'
#' @param iteration  Positive integer SLURM array task ID.
#' @param grid  Design grid tibble (default: DESIGN_GRID).
#' @param nc,nt  Sample sizes.
#' @param K_tt  Treated-to-treated neighbours for alt_tt.
#' @return Tibble with 20 × 4 = 80 rows (20 cells × 4 estimators).
sim_master_multi <- function(iteration,
                              grid  = DESIGN_GRID,
                              nc    = 500,
                              nt    = 100,
                              K_tt  = 2) {

  make_seed <- function(iter, overlap_label, error_type, sigma1_extra) {
    overlap_idx <- match(as.character(overlap_label), OVERLAP_LABELS)
    error_idx   <- if (as.character(error_type) == "homo") 1L else 2L
    sigma1_idx  <- if (sigma1_extra == 0) 1L else 2L
    iter +
      1000L  * (overlap_idx - 1L) +
      10000L * (error_idx   - 1L) +
      50000L * (sigma1_idx  - 1L)
  }

  results <- vector("list", nrow(grid))

  for (j in seq_len(nrow(grid))) {
    row  <- grid[j, ]
    seed <- make_seed(iteration,
                      row$overlap_label, row$error_type, row$sigma1_extra)
    prop <- PROP_NC_UNIF[as.character(row$overlap_label)]

    start <- Sys.time()
    results[[j]] <- tryCatch(
      one_iter(
        i             = iteration,
        overlap_label = as.character(row$overlap_label),
        error_type    = as.character(row$error_type),
        sigma1_extra  = row$sigma1_extra,
        prop_nc_unif  = prop,
        nc            = nc,
        nt            = nt,
        seed_addition = seed,
        K_tt          = K_tt
      ),
      error = function(e) {
        warning(sprintf(
          "one_iter failed: iter=%d overlap=%s error=%s sigma1=%.1f — %s",
          iteration, row$overlap_label, row$error_type,
          row$sigma1_extra, e$message
        ), call. = FALSE)
        tibble(
          runID            = iteration,
          inference_method = NA_character_,
          att_est          = NA_real_, SE = NA_real_,
          CI_lower         = NA_real_, CI_upper = NA_real_,
          V_E = NA_real_, V_P = NA_real_, SATT = NA_real_,
          overlap_label = as.character(row$overlap_label),
          error_type    = as.character(row$error_type),
          sigma1_extra  = row$sigma1_extra,
          common_label  = as.character(row$common_label),
          status        = "Iteration Failed"
        )
      }
    )
    results[[j]]$time_secs <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  }

  bind_rows(results)
}

# ─────────────────────────────────────────────────────────────────────────────
# Path helpers
# ─────────────────────────────────────────────────────────────────────────────

get_sim_paths <- function(output_name = "sims-variance-multi") {
  base <- here::here("data/outputs", output_name)
  dir.create(file.path(base, "individual"), showWarnings = FALSE, recursive = TRUE)
  list(
    base           = base,
    individual_dir = file.path(base, "individual"),
    combined_rds   = file.path(base, "combined_results.rds"),
    combined_csv   = file.path(base, "combined_results.csv")
  )
}
