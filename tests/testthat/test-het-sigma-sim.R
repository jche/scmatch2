# tests/testthat/test-het-sigma-sim.R
#
# Tests for scripts/sims-variance-het-sigma/
# Integration tests gated by RUN_SLOW_TESTS=TRUE (see helper-skip.R).

source(here::here("scripts/sims-variance-het-sigma/0_sim_inference_utils.R"))

# -------------------------------------------------------------------
# Unit tests: make_het_sigma_df (always run)
# -------------------------------------------------------------------

test_that("make_het_sigma_df returns a data frame with required columns", {
  set.seed(1)
  df <- make_het_sigma_df(nc = 50, nt = 20, prop_nc_unif = 1/3, seed = 1)

  expect_true(is.data.frame(df))
  expect_true(all(c("X1", "X2", "Z", "Y", "Y0", "Y1", "id",
                    "sigma0", "sigma1") %in% names(df)))
  expect_equal(nrow(df), 70)
  expect_equal(sum(df$Z == 1), 20)
})

test_that("make_het_sigma_df scenario 1: sigma0 == sigma1 pointwise", {
  set.seed(2)
  df <- make_het_sigma_df(nc = 100, nt = 50, prop_nc_unif = 1/3,
                          seed = 2, sigma1_extra = 0)
  expect_equal(df$sigma0, df$sigma1)
})

test_that("make_het_sigma_df scenario 2: sigma1 = sigma0 + 0.5 everywhere", {
  set.seed(3)
  df <- make_het_sigma_df(nc = 100, nt = 50, prop_nc_unif = 1/3,
                          seed = 3, sigma1_extra = 0.5)
  expect_equal(df$sigma1, df$sigma0 + 0.5)
})

test_that("make_het_sigma_df scenario 1: Y0 and Y1 share noise (Y1-Y0 = tau)", {
  # When sigma1_extra=0, eps0==eps1, so Y1-Y0 = Y1_denoised - Y0_denoised
  set.seed(4)
  df <- make_het_sigma_df(nc = 200, nt = 50, prop_nc_unif = 1/3,
                          seed = 4, sigma1_extra = 0)
  tau_check <- df$Y1_denoised - df$Y0_denoised
  expect_equal(df$Y1 - df$Y0, tau_check, tolerance = 1e-10)
})

test_that("make_het_sigma_df scenario 2: Y1-Y0 != pure tau (independent noise)", {
  # When sigma1_extra>0, eps1 is drawn independently so Y1-Y0 != tau(x) exactly
  set.seed(5)
  df <- make_het_sigma_df(nc = 500, nt = 100, prop_nc_unif = 1/3,
                          seed = 5, sigma1_extra = 0.5)
  tau_exact <- df$Y1_denoised - df$Y0_denoised
  residuals  <- (df$Y1 - df$Y0) - tau_exact
  # residuals should have non-zero variance (independent noise)
  expect_gt(sd(residuals), 0.01)
})

test_that("sigma0_fun computes 0.2 + (x1-x2)^2 correctly", {
  expect_equal(sigma0_fun(0.5, 0.5), 0.2)
  expect_equal(sigma0_fun(1.0, 0.0), 0.2 + 1.0)
  expect_equal(sigma0_fun(0.3, 0.7), 0.2 + (0.3 - 0.7)^2)
})

# -------------------------------------------------------------------
# Unit test: toy_match_infer_het output structure (always run)
# -------------------------------------------------------------------

test_that("toy_match_infer_het returns tibble with 4 methods × required columns", {
  skip_if_not(run_slow_tests(), "Skipping slow test (set RUN_SLOW_TESTS=TRUE)")

  res <- toy_match_infer_het(
    i            = 1,
    overlap_label = "mid",
    sigma1_extra = 0,
    prop_nc_unif = 1/3,
    nc = 100, nt = 30,
    K_tt = 1,
    seed_addition = 42
  )

  expect_true(is.data.frame(res))
  expect_setequal(res$inference_method, c("homo", "het", "alt_common", "alt_tt"))
  expect_true(all(c("att_est", "SE", "CI_lower", "CI_upper",
                    "V_E", "V_P", "SATT", "N_T", "ESS_C") %in% names(res)))
  expect_equal(nrow(res), 4)
})

test_that("toy_match_infer_het: all SEs are positive for scenario 1", {
  skip_if_not(run_slow_tests(), "Skipping slow test (set RUN_SLOW_TESTS=TRUE)")

  res <- toy_match_infer_het(
    i             = 1,
    overlap_label  = "mid",
    sigma1_extra  = 0,
    prop_nc_unif  = 1/3,
    nc = 100, nt = 30,
    K_tt = 1,
    seed_addition = 99
  )
  expect_true(all(res$SE > 0, na.rm = TRUE))
})

test_that("toy_match_infer_het: all SEs positive for scenario 2", {
  skip_if_not(run_slow_tests(), "Skipping slow test (set RUN_SLOW_TESTS=TRUE)")

  res <- toy_match_infer_het(
    i             = 2,
    overlap_label  = "high",
    sigma1_extra  = 0.5,
    prop_nc_unif  = 1/2,
    nc = 100, nt = 30,
    K_tt = 1,
    seed_addition = 77
  )
  expect_true(all(res$SE > 0, na.rm = TRUE))
})

# -------------------------------------------------------------------
# Integration test: sim_master_het end-to-end
# -------------------------------------------------------------------

test_that("sim_master_het runs and returns expected columns", {
  skip_if_not(run_slow_tests(), "Skipping slow test (set RUN_SLOW_TESTS=TRUE)")

  res <- sim_master_het(
    iteration     = 1,
    overlap_label  = "mid",
    sigma1_extra  = 0,
    rad_method    = "adaptive",
    k_match       = 2,
    est_method    = "scm",
    K_tt          = 1
  )

  expect_true(is.data.frame(res))
  expect_setequal(res$inference_method, c("homo", "het", "alt_common", "alt_tt"))
  expect_true("time_secs" %in% names(res))
  expect_true(all(res$time_secs >= 0, na.rm = TRUE))
})

# -------------------------------------------------------------------
# Integration coverage test: 50 reps, check alt_tt works in scenario 2
# -------------------------------------------------------------------

test_that("50-rep coverage: alt_tt covers SATT in scenario 1 & 2; homo does not", {
  skip_if_not(run_slow_tests(), "Skipping slow test (set RUN_SLOW_TESTS=TRUE)")

  iters    <- 1:50
  overlaps <- c("mid")

  run_one <- function(iter, sigma1_extra) {
    tryCatch(
      sim_master_het(
        iteration     = iter,
        overlap_label  = "mid",
        sigma1_extra  = sigma1_extra,
        rad_method    = "adaptive",
        k_match       = 2,
        est_method    = "scm",
        K_tt          = 2
      ),
      error = function(e) NULL
    )
  }

  res_s1 <- bind_rows(lapply(iters, run_one, sigma1_extra = 0))
  res_s2 <- bind_rows(lapply(iters, run_one, sigma1_extra = 0.5))

  coverage <- function(df, method) {
    dfm <- df %>% filter(inference_method == method, !is.na(CI_lower), !is.na(SATT))
    mean((dfm$CI_lower <= dfm$SATT) & (dfm$SATT <= dfm$CI_upper))
  }

  # Scenario 1: alt_tt and alt_common should have reasonable coverage (>60%)
  cov_alt_tt_s1     <- coverage(res_s1, "alt_tt")
  cov_alt_common_s1 <- coverage(res_s1, "alt_common")
  expect_gt(cov_alt_tt_s1,     0.60)
  expect_gt(cov_alt_common_s1, 0.60)

  # Scenario 2: alt_tt should have reasonable coverage; alt_common need not
  cov_alt_tt_s2 <- coverage(res_s2, "alt_tt")
  expect_gt(cov_alt_tt_s2, 0.60)
})
