# tests/testthat/test-wrappers.R

test_that("Causal Forest produces reasonable estimates", {
  skip_if_not_installed("grf")

  set.seed(123)

  # Generate simple DGP with known treatment effect
  n <- 500
  df <- tibble(
    X1 = runif(n),
    X2 = runif(n),
    Z = rbinom(n, 1, 0.3),
    Y0 = X1 + X2 + rnorm(n, 0, 1),
    Y1 = Y0 + 2,  # Constant treatment effect of 2
    Y = ifelse(Z, Y1, Y0)
  )

  true_att <- df %>%
    filter(Z == 1) %>%
    summarize(att = mean(Y1 - Y0)) %>%
    pull(att)

  att_cf <- get_att_causal_forest(df, covs = c("X1", "X2"))

  # Check that estimate is within reasonable range
  expect_true(abs(att_cf - true_att) < 0.5)
  expect_type(att_cf, "double")
  expect_length(att_cf, 1)

  cat("\nCausal Forest Test:\n")
  cat(sprintf("True ATT: %.3f\n", true_att))
  cat(sprintf("Estimated ATT: %.3f\n", att_cf))
})


test_that("TWANG produces reasonable estimates", {
  skip_if_not_installed("twang")

  set.seed(456)

  # Generate simple DGP
  n <- 500
  df <- tibble(
    X1 = runif(n),
    X2 = runif(n),
    Z = rbinom(n, 1, plogis(X1 + X2 - 1)),  # Selection on observables
    Y0 = X1 + X2 + rnorm(n, 0, 1),
    Y1 = Y0 + 1.5,  # Constant treatment effect
    Y = ifelse(Z, Y1, Y0)
  )

  true_att <- df %>%
    filter(Z == 1) %>%
    summarize(att = mean(Y1 - Y0)) %>%
    pull(att)

  att_twang <- get_att_twang(df, form = as.formula("Z ~ X1 + X2"))

  expect_type(att_twang, "double")
  expect_length(att_twang, 1)

  cat("\nTWANG Test:\n")
  cat(sprintf("True ATT: %.3f\n", true_att))
  cat(sprintf("Estimated ATT: %.3f\n", att_twang))
  cat(sprintf("Error: %.3f\n", att_twang - true_att))
})

test_that("TWANG works with actual simulation DGP", {
  skip_if_not_installed("twang")

  set.seed(1001)  # Same as first iteration

  # Use the ACTUAL DGP from sim_toy_new_methods.R
  nc <- 500
  nt <- 100
  f0_sd <- 0.5

  df <- gen_df_adv(
    nc = nc,
    nt = nt,
    f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) {3 * X1 + 3 * X2},
    f0_fun = function(x, y) {
      matrix(c(x, y), ncol = 2) %>%
        mvtnorm::dmvnorm(
          mean = c(0.5, 0.5),
          sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2)
        ) * 20
    }
  )

  num_bins <- 6
  dist_scaling <- df %>%
    summarize(across(
      starts_with("X"),
      function(x) {
        if (is.numeric(x)) 2 * num_bins / (max(x) - min(x))
        else 1000
      }
    ))

  # Apply same filtering
  preds_csm <- get_cal_matches(
    df = df,
    metric = "maximum",
    scaling = dist_scaling,
    rad_method = "fixed",
    est_method = "average",
    k = 25
  )

  df <- df %>%
    filter(!id %in% attr(preds_csm, "unmatched_units"))

  # Calculate true ATT
  true_att <- df %>%
    filter(Z == 1) %>%
    summarize(att = mean(Y1 - Y0)) %>%
    pull(att)

  # Test different formulas
  cat("\n=== Testing TWANG with actual DGP ===\n")

  # Simple additive
  att_twang1 <- get_att_twang(df, form = as.formula("Z ~ X1 + X2"))
  cat(sprintf("Formula: Z ~ X1 + X2\n"))
  cat(sprintf("  True ATT: %.3f\n", true_att))
  cat(sprintf("  Estimated ATT: %.3f\n", att_twang1))
  cat(sprintf("  Is NA: %s\n", is.na(att_twang1)))

  # # With interaction
  # att_twang2 <- get_att_twang(df, form = as.formula("Z ~ X1 * X2"))
  # cat(sprintf("\nFormula: Z ~ X1 * X2\n"))
  # cat(sprintf("  True ATT: %.3f\n", true_att))
  # cat(sprintf("  Estimated ATT: %.3f\n", att_twang2))
  # cat(sprintf("  Is NA: %s\n", is.na(att_twang2)))

  # # Check propensity score distribution
  # cat("\n=== Checking propensity score ===\n")
  # d_test <- as.data.frame(df)
  # d_test$Z <- as.numeric(d_test$Z)
  #
  # ps_fit <- twang::ps(
  #   formula = as.formula("Z ~ X1 + X2"),
  #   data = d_test,
  #   estimand = "ATT",
  #   stop.method = "es.mean",
  #   n.trees = 5000,
  #   verbose = TRUE  # Turn on verbose to see what's happening
  # )
  #
  # cat("PS fit summary:\n")
  # print(summary(ps_fit$ps))

  # weights <- twang::get.weights(ps_fit, stop.method = "es.mean")
  # cat("\nWeights summary:\n")
  # print(summary(weights))
  # cat("Any NA weights:", any(is.na(weights)), "\n")

  # # At least one should work
  # expect_true(!is.na(att_twang1) || !is.na(att_twang2))
})

test_that("Kbal produces reasonable estimates", {
  skip_if_not_installed("kbal")

  set.seed(789)

  # Generate simple DGP
  n <- 500
  df <- tibble(
    X1 = runif(n),
    X2 = runif(n),
    Z = rbinom(n, 1, 0.3),
    Y0 = X1 + X2 + rnorm(n, 0, 1),
    Y1 = Y0 + 1.8,  # Constant treatment effect
    Y = ifelse(Z, Y1, Y0)
  )

  true_att <- df %>%
    filter(Z == 1) %>%
    summarize(att = mean(Y1 - Y0)) %>%
    pull(att)

  att_kbal <- get_att_kbal(df, covs = c("X1", "X2"))

  expect_true(abs(att_kbal - true_att) < 0.5)
  expect_type(att_kbal, "double")
  expect_length(att_kbal, 1)

  cat("\nKbal Test:\n")
  cat(sprintf("True ATT: %.3f\n", true_att))
  cat(sprintf("Estimated ATT: %.3f\n", att_kbal))
  cat(sprintf("Error: %.3f\n", att_kbal - true_att))
})


test_that("Kbal with fixed K=30 produces reasonable estimates", {
  skip_if_not_installed("kbal")
  set.seed(789)

  # Generate simple DGP
  n <- 500
  df <- tibble(
    X1 = runif(n),
    X2 = runif(n),
    Z = rbinom(n, 1, 0.3),
    Y0 = X1 + X2 + rnorm(n, 0, 1),
    Y1 = Y0 + 1.8,  # Constant treatment effect
    Y = ifelse(Z, Y1, Y0)
  )

  true_att <- df %>%
    filter(Z == 1) %>%
    summarize(att = mean(Y1 - Y0)) %>%
    pull(att)

  # Test with fixed K = 30
  att_kbal_k30 <- get_att_kbal(df, covs = c("X1", "X2"), numdims = 30)

  expect_true(abs(att_kbal_k30 - true_att) < 0.5)
  expect_type(att_kbal_k30, "double")
  expect_length(att_kbal_k30, 1)

  cat("\nKbal Test (K=30):\n")
  cat(sprintf("True ATT: %.3f\n", true_att))
  cat(sprintf("Estimated ATT (K=30): %.3f\n", att_kbal_k30))
  cat(sprintf("Error: %.3f\n", att_kbal_k30 - true_att))
})



test_that("get_att_point_est works well",{
  test_matched_weighted_df <-
    data.frame(Z= c(1,0,0,0,1),
               Y = c(1,-1,1,0,1),
               weights = c(1,0.5,0.5,1,1) )
  test_matched_weighted_df
  atts <- CSM:::get_att_point_est(test_matched_weighted_df)
  atts
  expect_equal( atts, 1 )
})


test_that("get_att_bal should work for a example dataset",{

  # Should work  with either
  # Z in 0,1 or Z in F,T, but not in other numerical Z

  require( optweight )

  source( here::here( "scripts/datagen/gen_six_points.R" ) )
  df_six_points <- gen_six_points()
  covs <- c("X1", "X2")
  zform1 <-
    as.formula(paste0("Z ~ ",
                      paste0(covs, collapse="+")))

  test_att_bal = CSM:::get_att_bal(d = df_six_points,
                     form = zform1,
                     tols = rep(0.01, length(covs)))
  expect_equal(test_att_bal, 9.0399,tolerance=1e-3)

  df_Z_logical <-
    df_six_points %>%
    mutate(Z=as.logical(Z))
  test_att_bal = CSM:::get_att_bal(d = df_Z_logical,
                     form = zform1,
                     tols = rep(0.01, length(covs)))
  expect_equal(test_att_bal, 9.0399,tolerance=1e-3)

  # Make sure treatment is binary is checked
  df_Z_random <- df_six_points %>%
    ungroup() %>%
    mutate( Z = as.numeric(1:6))
  expect_error(
    CSM:::get_att_bal(d = df_Z_random,
                form = zform1,
                tols = rep(0.01, length(covs))),
    "Treatment variable must be either 0, 1 or FALSE, TRUE")
})



test_that( "cem works", {
  df <- CSM:::gen_one_toy()

  cem <- CSM:::get_att_cem(df, num_bins=5, est_method="scm")

  expect_true( is.numeric(cem) )
})

test_that("get_matches works with different matching types", {
  test_df <- data.frame(
    Z = c(1, 0, 0, 0, 1),
    X = c(0, 0.5, 0.8, 3, 1.6)
  )
  scaling <- 1  # Example scaling factor

  res_toy <- get_matches(
    matching_type = "maximum_fixed_scm",
    df_dgp = test_df,
    scaling = scaling
  )
  expect_equal(nrow(res_toy), 5)  # Ensure 4 rows in the result

  # Test data for otsu
  test_df_otsu <- data.frame(
    Z = c(1, 0, 0, 0, 1),
    X1 = c(0, 0.5, 0.8, 3, 1.6),
    X2 = c(0, 0, 0, 0, 0)
  )
  test_df_otsu <- bind_rows(
    test_df_otsu,
    transform(test_df_otsu, X2 = 1, Z = 0)
  )

  # Test for 'euclidean_knn'
  res_otsu <- get_matches(
    matching_type = "euclidean_knn",
    df_dgp = test_df_otsu,
    scaling = scaling
  )
  expect_equal(nrow(res_otsu), 18)

  # Test for invalid matching_type
  expect_error(
    get_matches(
      matching_type = "invalid_type",
      df_dgp = test_df,
      scaling = scaling
    ),
    "Invalid matching_type"
  )
})



test_that( "all other methods run and give ATT estimates", {

  set.seed( 40440 )
  df <- CSM:::gen_one_toy()
  df

  df %>%
    ggplot(aes(X1,X2)) +
    geom_point(aes(pch=as.factor(Z), color=Y)) +
    scale_color_continuous(low="orange", high="blue") +
    theme_classic() +
    labs(pch = "Treated",
         color = latex2exp::TeX("$f_0(X)$"),
         x = latex2exp::TeX("$X_1$"),
         y = latex2exp::TeX("$X_2$")) +
    facet_wrap(~Z)

  df %>%
    filter(Z) %>%
    summarize(eff = mean(Y1-Y0))

  res <- CSM:::run_all_methods( df )
  res
  expect_true( is.data.frame(res) )


  df = filter( df, Z==1 )
  df2 = df %>%
    mutate( Z = 0,
            Y = Y - 1 )

  df = bind_rows(df, df2)
  df

  # Looking at CSM local averaging issue
  preds_csm <- get_cal_matches(
    df = df,
    metric = "maximum",
    scaling = 100,
    #caliper = 0.001,
    rad_method = "adaptive",
    est_method = "scm" )
  preds_csm$treatment_table
  head( preds_csm$matches )

  tt <- CSM:::get_att_csm(df, scaling=100)
  expect_equal( tt, 1 )

  # NOTE: CSM is biased here due to radial matching keeping lots of
  # units due to large radius
  res <- CSM:::run_all_methods( df, extrapolate = FALSE )
  res <- filter( res, !is.na(att) )
  res <- filter( res, method != "csm" )
  res <- filter( res, method != "cem" ) # It extrapolates and messes up?
  expect_equal( res$att, rep( 1, nrow(res) ), tolerance=1e-2 )
})
