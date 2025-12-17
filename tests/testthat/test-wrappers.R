# tests/testthat/test-wrappers.R
source(here::here("scripts/wrappers.R"))

## ---------------------- FAST TESTS (always run) ----------------------

test_that("Causal Forest produces reasonable estimates (small, fast DGP)", {
  skip_if_not_installed("grf")

  set.seed(123)
  n <- 300
  df <- tibble::tibble(
    X1 = runif(n),
    X2 = runif(n),
    Z  = rbinom(n, 1, 0.3),
    Y0 = X1 + X2 + rnorm(n, 0, 1),
    Y1 = Y0 + 2,
    Y  = ifelse(Z, Y1, Y0)
  )

  true_att <- df |>
    dplyr::filter(Z == 1) |>
    dplyr::summarize(att = mean(Y1 - Y0)) |>
    dplyr::pull(att)

  att_cf <- get_att_causal_forest(df, covs = c("X1", "X2"))
  expect_true(abs(att_cf - true_att) < 0.6)
  expect_type(att_cf, "double")
  expect_length(att_cf, 1)
})

test_that("get_att_point_est works", {
  dat <- data.frame(Z = c(1,0,0,0,1),
                    Y = c(1,-1,1,0,1),
                    weights = c(1,0.5,0.5,1,1))
  atts <- CSM:::get_att_point_est(dat)
  expect_equal(atts, 1)
})


test_that("cem returns numeric", {
  df <- CSM:::gen_one_toy()
  cem <- get_att_cem(df, num_bins = 5, est_method = "scm")
  expect_true(is.numeric(cem))
})



test_that("get_matches works with different matching types", {
  test_df <- data.frame(Z = c(1,0,0,0,1),
                        X = c(0,0.5,0.8,3,1.6))
  res_max <- get_matches("maximum_fixed_scm",  test_df, scaling = 1)
  expect_true(nrow(res_max) >= 1)

  # euclidean_knn minimal check
  test_df2 <- data.frame(
    Z = c(1,0,0,0,1),
    X1 = c(0,0.5,0.8,3,1.6),
    X2 = c(0,0,0,0,0)
  )
  res_knn <- get_matches("euclidean_knn", test_df2, scaling = 1)
  expect_true(nrow(res_knn) >= 1)

  expect_error(
    get_matches("invalid_type", test_df, scaling = 1),
    "Invalid matching_type"
  )
})

## ---------------------- SLOW TESTS (skipped by default) ----------------------

test_that("tmle2 and aipw2 complete on toy (slow)", {
  skip_if_not(run_slow_tests(), "Set RUN_SLOW_TESTS=TRUE to enable slow tests")
  skip_if_not_installed("tmle")
  skip_if_not_installed("dbarts")
  skip_if_not_installed("AIPW")

  set.seed(42)
  df <- CSM:::gen_one_toy(nc = 120, nt = 40, f0_sd = 0.2)  # keep small
  covs <- c("X1","X2")

  # These call heavier SL libraries internally; keep dataset tiny
  tmle2 <- get_att_tmle(df, covs = covs,
                        Q.SL.library = c("SL.glm","tmle.SL.dbarts2","SL.glmnet"),
                        g.SL.library = c("SL.glm","tmle.SL.dbarts.k.5","SL.gam"))
  aipw2 <- get_att_aipw(df, covs = covs,
                        Q.SL.library = c("SL.glm","SL.glmnet","SL.randomForest","SL.xgboost"),
                        g.SL.library = c("SL.glm","SL.glmnet","SL.randomForest","SL.xgboost"))
  expect_true(is.numeric(tmle2) || is.numeric(aipw2))
})

test_that("twang returns numeric for hainmueller (slow)", {
  skip_if_not(run_slow_tests(), "Set RUN_SLOW_TESTS=TRUE to enable slow tests")
  skip_if_not_installed("twang")

  set.seed(123)
  df <- CSM:::gen_df_hain(nc = 200, nt = 40, sigma_e = "n100", outcome = "nl2", sigma_y = 1, ATE = 0)
  est <- get_att_twang(df, form = as.formula("Z ~ X1+X2+X3+X4+X5+X6"))
  expect_true(is.numeric(est) && length(est) == 1)
  expect_false(is.na(est))
})

test_that("TWANG synthetic DGP", {
  #skip_if_not(run_slow_tests(), "Set RUN_SLOW_TESTS=TRUE to enable slow tests")
  skip_if_not_installed("twang")

  set.seed(456)
  n <- 400
  df <- tibble::tibble(
    X1 = runif(n),
    X2 = runif(n),
    Z  = rbinom(n, 1, plogis(X1 + X2 - 1)),
    Y0 = X1 + X2 + rnorm(n, 0, 1),
    Y1 = Y0 + 1.5,
    Y  = ifelse(Z, Y1, Y0)
  )

  true_att <- df |>
    dplyr::filter(Z == 1) |>
    dplyr::summarize(att = mean(Y1 - Y0)) |>
    dplyr::pull(att)

  att_twang <- get_att_twang(df, form = as.formula("Z ~ X1 + X2"))
  expect_type(att_twang, "double")
  expect_length(att_twang, 1)
})

test_that("Kbal estimates", {
 # skip_if_not(run_slow_tests(), "Set RUN_SLOW_TESTS=TRUE to enable slow tests")
  skip_if_not_installed("kbal")

  set.seed(789)
  n <- 400
  df <- tibble::tibble(
    X1 = runif(n),
    X2 = runif(n),
    Z  = rbinom(n, 1, 0.3),
    Y0 = X1 + X2 + rnorm(n, 0, 1),
    Y1 = Y0 + 1.8,
    Y  = ifelse(Z, Y1, Y0)
  )
  att_kbal <- get_att_kbal(df, covs = c("X1","X2"), numdims = 20) # smaller numdims
  expect_type(att_kbal, "double")
  expect_length(att_kbal, 1)
})
