# scripts/figs/tune_bimodal_sigma_params.R
#
# Grid search over (A, h_sq, sigma_min) for the bimodal σ₀ design.
# Goal: find parameters that produce a much larger Cov_p(w_j, s_j²),
# so that the coverage gap between homo and het estimators is clearly visible.
#
# Current baseline: A=1, h_sq=0.08, sigma_min=0.2 → Cov_p ≈ +0.115

devtools::load_all()

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(here)
  library(mvtnorm)
})

source(here::here("scripts/sims-variance-het-sigma-bimodal/0_sim_inference_utils.R"))

# -------------------------------------------------------------------
# Parameterized σ₀ (same form as bimodal but freely tunable)
# -------------------------------------------------------------------

sigma0_tunable <- function(x1, x2, sigma_min, A, h_sq) {
  d1_sq <- (x1 - 0.25)^2 + (x2 - 0.25)^2
  d2_sq <- (x1 - 0.75)^2 + (x2 - 0.75)^2
  sigma_min + A * (exp(-d1_sq / h_sq) + exp(-d2_sq / h_sq))
}

# -------------------------------------------------------------------
# Single-dataset Cov_p + V_E estimates for one parameter set
# -------------------------------------------------------------------

eval_params_once <- function(A, h_sq, sigma_min, seed) {
  set.seed(seed)

  df_raw <- gen_df_adv_k(
    nc = 500, nt = 100, k = 2,
    f0_sd = 0,
    tx_effect_fun = function(X) {
      X <- as.matrix(X); 3 * X[, 1] + 3 * X[, 2]
    },
    f0_fun = function(X) {
      X <- as.matrix(X)
      mvtnorm::dmvnorm(X, mean  = c(0.5, 0.5),
                       sigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)) * 20
    },
    ctr_dist = 0.5, prop_nc_unif = 1/3
  )

  sig0 <- sigma0_tunable(df_raw$X1, df_raw$X2,
                         sigma_min = sigma_min, A = A, h_sq = h_sq)
  df <- df_raw %>%
    mutate(
      Y0 = Y0_denoised + rnorm(n(), 0, sig0),
      Y1 = Y1_denoised + rnorm(n(), 0, sig0),
      Y  = ifelse(Z, Y1, Y0),
      Z  = as.integer(Z)
    )

  scaling <- compute_toy_scaling(df)

  mtch <- tryCatch(
    get_cal_matches(data = df, Z ~ X1 + X2,
                    rad_method = "adaptive", scaling = scaling,
                    k = 2, warn = FALSE, est_method = "scm"),
    error = function(e) NULL
  )
  if (is.null(mtch)) return(NULL)

  # --- Cov_p(w_j, s_j²) ---
  co_wts <- result_table(mtch, return = "agg_co_units") %>%
    filter(Z == 0) %>% select(id, w_j = weights)
  df_ctrl <- df %>% filter(Z == 0) %>% mutate(id = as.character(id))
  mf  <- full_unit_table(mtch)
  sj  <- calculate_s_j_sq(mf, outcome = "Y", treatment = "Z")

  plot_df <- df_ctrl %>%
    left_join(co_wts, by = "id") %>%
    mutate(w_j = replace(w_j, is.na(w_j), 0)) %>%
    left_join(sj$s_j_sq %>% mutate(id = as.character(id)), by = "id") %>%
    filter(w_j > 1e-9, !is.na(s_j_sq))

  cov_p <- with(plot_df, cov(w_j, s_j_sq))

  # --- V_E: homo vs het ---
  V_E_homo <- tryCatch(
    estimate_ATT(mtch, variance_method = "pooled")$V_E,
    error = function(e) NA_real_
  )
  V_E_het <- tryCatch(
    estimate_ATT(mtch, variance_method = "pooled_het")$V_E,
    error = function(e) NA_real_
  )

  tibble(cov_p = cov_p, V_E_homo = V_E_homo, V_E_het = V_E_het)
}

# -------------------------------------------------------------------
# Grid
# -------------------------------------------------------------------

grid <- expand.grid(
  A         = c(1, 2, 3, 4, 5),
  h_sq      = c(0.04, 0.08, 0.12),
  sigma_min = c(0.05, 0.1, 0.2),
  stringsAsFactors = FALSE
)

n_reps  <- 5
base_seed <- 42

message("Running ", nrow(grid), " × ", n_reps, " = ", nrow(grid) * n_reps,
        " evaluations …")

results <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  A         <- grid$A[i]
  h_sq      <- grid$h_sq[i]
  sigma_min <- grid$sigma_min[i]

  reps <- lapply(seq_len(n_reps), function(r) {
    eval_params_once(A, h_sq, sigma_min, seed = base_seed + r * 100 + i)
  })
  reps <- bind_rows(Filter(Negate(is.null), reps))

  results[[i]] <- tibble(
    A         = A,
    h_sq      = h_sq,
    sigma_min = sigma_min,
    cov_p     = mean(reps$cov_p,    na.rm = TRUE),
    V_E_homo  = mean(reps$V_E_homo, na.rm = TRUE),
    V_E_het   = mean(reps$V_E_het,  na.rm = TRUE),
    VE_ratio  = mean(reps$V_E_het,  na.rm = TRUE) /
                mean(reps$V_E_homo, na.rm = TRUE),
    n_ok      = sum(!is.na(reps$cov_p))
  )
  message(sprintf("  A=%.0f  h_sq=%.2f  sigma_min=%.2f  →  Cov_p=%.3f  VE_ratio=%.3f",
                  A, h_sq, sigma_min,
                  results[[i]]$cov_p, results[[i]]$VE_ratio))
}

res_df <- bind_rows(results) %>%
  arrange(desc(cov_p))

cat("\n=== Results sorted by Cov_p (descending) ===\n")
print(as.data.frame(res_df), digits = 3)

# -------------------------------------------------------------------
# Highlight the "best" candidates (large Cov_p, large VE_ratio)
# -------------------------------------------------------------------

cat("\n=== Top-5 by VE_ratio ===\n")
print(as.data.frame(res_df %>% arrange(desc(VE_ratio)) %>% head(5)), digits = 3)

out_path <- here::here("scripts/figs/tune_bimodal_sigma_params.rds")
saveRDS(res_df, out_path)
message("\nSaved: ", out_path)
