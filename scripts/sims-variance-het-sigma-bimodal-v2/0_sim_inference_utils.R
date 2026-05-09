# scripts/sims-variance-het-sigma-bimodal-v2/0_sim_inference_utils.R
#
# Like sims-variance-het-sigma-bimodal/ (v1), but with stronger heteroskedasticity:
#
#   A=5, h_sq=0.08, sigma_min=0.2
#   σ₀(x) = 0.2 + 5·[exp(-‖x-(0.25,0.25)‖²/0.08) + exp(-‖x-(0.75,0.75)‖²/0.08)]
#
# At the treated cluster centres: σ₀ ≈ 5.2  (σ₀² ≈ 27)
# At the control cluster centres: σ₀ ≈ 0.2  (σ₀² ≈ 0.04)
#
# This amplifies Cov_p(w_j, s_j²) from ≈0.115 (v1) to ≈4.65,
# making the V_E_het / V_E_homo ratio ≈ 1.20 (vs ≈ 1.04 in v1).
#
# Scenarios:
#   sigma1_extra=0   → σ₁=σ₀  (homo should mildly undercover; het should be nominal)
#   sigma1_extra=0.5 → σ₁=σ₀+0.5  (only alt_tt should work)
#
# Parameter choice rationale: grid search over A∈{1..5}, h_sq∈{0.04,0.08,0.12},
# sigma_min∈{0.05,0.1,0.2} run in scripts/figs/tune_bimodal_sigma_params.R.

devtools::load_all()

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(here)
  library(mvtnorm)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (is.atomic(x) && all(is.na(x)))) y else x
}

get_sim_paths <- function(output_name = "sims-variance-het-sigma-bimodal-v2") {
  output_dir_base <- here::here(file.path("data/outputs", output_name))
  dir.create(output_dir_base, showWarnings = FALSE, recursive = TRUE)
  list(
    output_dir_base = output_dir_base,
    individual_dir  = file.path(output_dir_base, "individual"),
    summary_csv     = file.path(output_dir_base, "summary_table.csv"),
    combined_csv    = file.path(output_dir_base, "combined_results.csv")
  )
}

# -------------------------------------------------------------------
# σ₀(x): bimodal Gaussian bumps, v2 parameters (A=5, h_sq=0.08, sigma_min=0.2)
# -------------------------------------------------------------------

SIGMA0_V2_PARAMS <- list(sigma_min = 0.2, A = 5.0, h_sq = 0.08)

sigma0_bimodal_v2 <- function(x1, x2,
                               sigma_min = SIGMA0_V2_PARAMS$sigma_min,
                               A         = SIGMA0_V2_PARAMS$A,
                               h_sq      = SIGMA0_V2_PARAMS$h_sq) {
  d1_sq <- (x1 - 0.25)^2 + (x2 - 0.25)^2
  d2_sq <- (x1 - 0.75)^2 + (x2 - 0.75)^2
  sigma_min + A * (exp(-d1_sq / h_sq) + exp(-d2_sq / h_sq))
}

# -------------------------------------------------------------------
# Data generation
# -------------------------------------------------------------------

make_bimodal_sigma_v2_df <- function(
    nc = 500, nt = 100,
    prop_nc_unif = 1/3,
    ctr_dist     = 0.5,
    seed         = NULL,
    sigma1_extra = 0
) {
  if (!is.null(seed)) set.seed(seed)

  f0_fun_mat <- function(X) {
    X <- as.matrix(X)
    mvtnorm::dmvnorm(X, mean  = c(0.5, 0.5),
                     sigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)) * 20
  }
  tx_fun_mat <- function(X) {
    X <- as.matrix(X)
    3 * X[, 1] + 3 * X[, 2]
  }

  df_raw <- gen_df_adv_k(
    nc = nc, nt = nt, k = 2,
    f0_sd         = 0,
    tx_effect_fun = tx_fun_mat,
    f0_fun        = f0_fun_mat,
    ctr_dist      = ctr_dist,
    prop_nc_unif  = prop_nc_unif
  )

  sig0 <- sigma0_bimodal_v2(df_raw$X1, df_raw$X2)
  sig1 <- sig0 + sigma1_extra

  eps0 <- rnorm(nrow(df_raw), 0, sig0)
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

# -------------------------------------------------------------------
# Scaling
# -------------------------------------------------------------------

compute_toy_scaling <- function(df, nbins = 6) {
  df %>%
    summarise(across(
      starts_with("X"),
      function(x) if (is.numeric(x)) (2 * nbins) / (max(x) - min(x)) else 1000
    ))
}

# -------------------------------------------------------------------
# One iteration: match + four inference methods
# -------------------------------------------------------------------

toy_match_infer_bimodal_v2 <- function(
    i,
    overlap_label,
    sigma1_extra  = 0,
    prop_nc_unif,
    nc            = 500,
    nt            = 100,
    ctr_dist      = 0.5,
    nbins         = 6,
    rad_method    = "adaptive",
    k_match       = 2,
    est_method    = "scm",
    K_tt          = 2,
    seed_addition = 11,
    grid_id       = 1,
    verbose       = FALSE
) {
  set.seed(seed_addition)

  df <- make_bimodal_sigma_v2_df(
    nc = nc, nt = nt,
    prop_nc_unif = prop_nc_unif,
    ctr_dist     = ctr_dist,
    seed         = seed_addition,
    sigma1_extra = sigma1_extra
  )

  scaling <- compute_toy_scaling(df, nbins = nbins)

  mtch <- tryCatch(
    get_cal_matches(
      data       = df,
      Z ~ X1 + X2,
      rad_method = rad_method,
      scaling    = scaling,
      k          = k_match,
      warn       = FALSE,
      est_method = est_method
    ),
    error = function(e) {
      if (verbose) message("Matching failed iter ", i, ": ", e$message)
      NULL
    }
  )

  att_est   <- if (!is.null(mtch)) tryCatch(get_att_point_est(mtch), error = function(e) NA_real_) else NA_real_
  true_satt <- df %>% filter(Z == 1) %>% summarise(att = mean(Y1 - Y0)) %>% pull(att)

  bias_val <- NA_real_
  if (!is.null(mtch)) {
    full_units <- tryCatch(full_unit_table(mtch), error = function(e) NULL)
    if (!is.null(full_units) && all(c("Y0", "weights", "Z") %in% names(full_units))) {
      tmp <- full_units %>%
        group_by(Z) %>%
        summarise(mn = stats::weighted.mean(Y0, w = weights, na.rm = TRUE), .groups = "drop")
      if (nrow(tmp) == 2 && all(c(0L, 1L) %in% tmp$Z))
        bias_val <- tmp$mn[tmp$Z == 1] - tmp$mn[tmp$Z == 0]
    }
  }

  make_row <- function(method_name, SE, V_E, V_P) {
    tibble(
      runID = i, inference_method = method_name,
      att_est = att_est, SE = SE,
      CI_lower = att_est - 1.96 * SE, CI_upper = att_est + 1.96 * SE,
      V_E = V_E, V_P = V_P, SATT = true_satt, bias = bias_val
    )
  }
  make_na_row <- function(method_name) {
    tibble(
      runID = i, inference_method = method_name,
      att_est = att_est, SE = NA_real_,
      CI_lower = NA_real_, CI_upper = NA_real_,
      V_E = NA_real_, V_P = NA_real_, SATT = true_satt, bias = bias_val
    )
  }

  methods_all <- c("homo", "het", "alt_common", "alt_tt")

  if (is.null(mtch)) {
    return(bind_rows(lapply(methods_all, make_na_row)) %>%
             mutate(status = "Matching Failed"))
  }

  matches_full <- tryCatch(full_unit_table(mtch), error = function(e) NULL)

  row_homo <- tryCatch({
    res1 <- get_ATT_estimate(mtch, variance_method = "pooled")
    make_row("homo", res1$SE, res1$V_E, res1$V_P)
  }, error = function(e) { if (verbose) message("homo: ", e$message); make_na_row("homo") })

  row_het <- tryCatch({
    res2 <- get_ATT_estimate(mtch, variance_method = "pooled_het")
    make_row("het", res2$SE, res2$V_E, res2$V_P)
  }, error = function(e) { if (verbose) message("het: ", e$message); make_na_row("het") })

  V_P_ref <- row_homo$V_P %||% NA_real_

  row_alt_common <- tryCatch({
    if (is.null(matches_full)) stop("no matches_full")
    alt <- get_measurement_error_variance_alt(matches_full, use_common_variance = TRUE)
    SE_c <- sqrt(pmax(0, alt$V_E_alt + V_P_ref))
    make_row("alt_common", SE_c, alt$V_E_alt, V_P_ref)
  }, error = function(e) { if (verbose) message("alt_common: ", e$message); make_na_row("alt_common") })

  row_alt_tt <- tryCatch({
    if (is.null(matches_full)) stop("no matches_full")
    alt_tt <- get_measurement_error_variance_alt(
      matches_table = matches_full, df = df,
      use_common_variance = FALSE, K = K_tt,
      covs = c("X1", "X2"), scaling = scaling
    )
    SE_t <- sqrt(pmax(0, alt_tt$V_E_alt + V_P_ref))
    make_row("alt_tt", SE_t, alt_tt$V_E_alt, V_P_ref)
  }, error = function(e) { if (verbose) message("alt_tt: ", e$message); make_na_row("alt_tt") })

  rs <- bind_rows(row_homo, row_het, row_alt_common, row_alt_tt)

  sz <- tryCatch(
    calc_N_T_N_C(full_unit_table(mtch)),
    error = function(e) list(N_T = NA_real_, N_C_tilde = NA_real_)
  )

  rs %>%
    mutate(
      N_T          = sz$N_T,
      ESS_C        = sz$N_C_tilde,
      seed         = seed_addition,
      gridID       = grid_id,
      deg_overlap  = overlap_label,
      sigma1_extra = sigma1_extra,
      prop_nc_unif = prop_nc_unif,
      N = nc + nt, nc = nc, nt = nt,
      ctr_dist   = ctr_dist,
      nbins      = nbins,
      rad_method = rad_method,
      k_match    = k_match,
      est_method = est_method,
      status     = "Success"
    )
}

# -------------------------------------------------------------------
# sim_master_bimodal_v2
# -------------------------------------------------------------------

sim_master_bimodal_v2 <- function(
    iteration, N = 600, overlap_label,
    sigma1_extra = 0,
    rad_method = "adaptive", k_match = 2, est_method = "scm",
    k_dim = 2, grid_id = 1, K_tt = 2, ...
) {
  prop_nc_unif_values <- c(
    very_high = 2/3, high = 1/2, mid = 1/3, low = 1/5, very_low = 1/10
  )
  prop_nc_unif <- prop_nc_unif_values[[overlap_label]]
  if (is.null(prop_nc_unif)) stop("Unknown overlap_label: ", overlap_label)

  nt <- 100; nc <- 500

  overlap_index <- match(overlap_label,
                         c("very_low","low","mid","high","very_high"))
  if (is.na(overlap_index)) overlap_index <- 1L
  sigma1_index <- if (sigma1_extra == 0) 1L else 2L

  seed <- iteration + 1000L * (overlap_index - 1L) + 50000L * (sigma1_index - 1L)

  start <- Sys.time()
  result <- tryCatch(
    toy_match_infer_bimodal_v2(
      i = iteration, overlap_label = overlap_label,
      sigma1_extra = sigma1_extra,
      prop_nc_unif = prop_nc_unif,
      nc = nc, nt = nt, ctr_dist = 0.5, nbins = 6,
      rad_method = rad_method, k_match = k_match, est_method = est_method,
      K_tt = K_tt, seed_addition = seed, grid_id = grid_id
    ),
    error = function(e) {
      warning("sim_master_bimodal_v2 failed iter ", iteration,
              " overlap ", overlap_label, ": ", e$message, call. = FALSE)
      tibble(
        runID = iteration, gridID = grid_id,
        deg_overlap = overlap_label, sigma1_extra = sigma1_extra,
        prop_nc_unif = prop_nc_unif, N = nc + nt, nc = nc, nt = nt,
        att_est = NA_real_, SE = NA_real_, status = "Iteration Failed"
      )
    }
  )
  result$time_secs <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  result
}
