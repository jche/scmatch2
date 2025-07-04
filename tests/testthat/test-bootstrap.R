# source("R/bootstrap.R")

test_that("Bayesian bootstrap standard errors match theoretical expectations for Gaussian data", {
  set.seed(42)
  mean_X <- 0
  sd_X <- 2
  n <- 500
  I <- 100
  se_gaussian_expected <- sd_X / sqrt(n)
  tol <- 0.01

  # Test Bayesian method
  bootstrap_ci_bayesian <- make_bootstrap_ci("Bayesian")
  se_boot_bayesian <- covered_bayesian <- numeric(I)

  for (i in 1:I) {
    # Generate data
    X <- rnorm(n, mean = mean_X, sd = sd_X)
    # Calculate bootstrap SEs using Bayesian method
    results <- bootstrap_ci_bayesian(
      resids = X - mean(X),
      mean_est = mean(X),
      B = 250,
      seed_addition = i
    )
    se_boot_bayesian[i] <- results$sd
    covered_bayesian[i] <- (results$ci_lower < mean_X) & (mean_X < results$ci_upper)
  }

  # Assertions for Bayesian method
  expect_lt(
    abs(mean(se_boot_bayesian) - se_gaussian_expected),
    tol,
    label = "Bootstrap SE should approximate theoretical SE for Bayesian"
  )

  expect_true(
    all(se_boot_bayesian > 0),
    label = "All standard errors should be positive for Bayesian"
  )

  coverage_bayesian <- mean(covered_bayesian)
  expect_true(
    coverage_bayesian >= 0.9 && coverage_bayesian <= 1,
    label = "CI coverage should be between 90% and 100% for Bayesian"
  )
})

test_that("Wild bootstrap standard errors match theoretical expectations for Gaussian data", {
  set.seed(42)
  mean_X <- 0
  sd_X <- 2
  n <- 500
  I <- 100
  se_gaussian_expected <- sd_X / sqrt(n)
  tol <- 0.01

  # Test wild method
  bootstrap_ci_wild <- make_bootstrap_ci("wild")
  se_boot_wild <- covered_wild <- numeric(I)

  for (i in 1:I) {
    # Generate data
    X <- rnorm(n, mean = mean_X, sd = sd_X)
    # Calculate bootstrap SEs using wild method
    results <- bootstrap_ci_wild(
      resids = X - mean(X),
      mean_est = mean(X),
      B = 250,
      seed_addition = i
    )
    se_boot_wild[i] <- results$sd
    covered_wild[i] <- (results$ci_lower < mean_X) & (mean_X < results$ci_upper)
  }

  # Assertions for wild method
  expect_lt(
    abs(mean(se_boot_wild) - se_gaussian_expected),
    tol,
    label = "Bootstrap SE should approximate theoretical SE for wild"
  )

  expect_true(
    all(se_boot_wild > 0),
    label = "All standard errors should be positive for wild"
  )

  coverage_wild <- mean(covered_wild)
  expect_true(
    coverage_wild >= 0.9 && coverage_wild <= 1,
    label = "CI coverage should be between 90% and 100% for wild"
  )
})

test_that("Naive bootstrap standard errors match theoretical expectations for Gaussian data", {
  set.seed(42)
  mean_X <- 0
  sd_X <- 2
  n <- 500
  I <- 100
  se_gaussian_expected <- sd_X / sqrt(n)
  tol <- 0.01

  # Test naive method
  bootstrap_ci_naive <- make_bootstrap_ci("naive")
  se_boot_naive <- covered_naive <- numeric(I)

  for (i in 1:I) {
    # Generate data
    X <- rnorm(n, mean = mean_X, sd = sd_X)
    # Calculate bootstrap SEs using naive method
    results <- bootstrap_ci_naive(
      resids = X - mean(X),
      mean_est = mean(X),
      B = 250,
      seed_addition = i
    )
    se_boot_naive[i] <- results$sd
    covered_naive[i] <- (results$ci_lower < mean_X) & (mean_X < results$ci_upper)
  }

  # Assertions for naive method
  expect_lt(
    abs(mean(se_boot_naive) - se_gaussian_expected),
    tol,
    label = "Bootstrap SE should approximate theoretical SE for naive"
  )

  expect_true(
    all(se_boot_naive > 0),
    label = "All standard errors should be positive for naive"
  )

  coverage_naive <- mean(covered_naive)
  expect_true(
    coverage_naive >= 0.9 && coverage_naive <= 1,
    label = "CI coverage should be between 90% and 100% for naive"
  )
})


test_that("Moving Block Bootstrap standard errors match theoretical expectations for Gaussian data", {
  set.seed(42)
  mean_X <- 0
  sd_X <- 2
  n <- 500
  I <- 100
  se_gaussian_expected <- sd_X / sqrt(n)
  tol <- 0.05

  test_bootstrap_method <- function(method, mean_X, sd_X, n, I, se_gaussian_expected, tol) {
    # Create bootstrap CI calculator for the given method
    bootstrap_ci <- make_bootstrap_ci(method, use_moving_block=T)

    # Run simulation
    se_boot <- covered <- numeric(I)
    for (i in 1:I) {
      # Generate data
      X <- rnorm(n, mean = mean_X, sd = sd_X)

      # Calculate bootstrap SEs using the selected method
      results <- bootstrap_ci(
        resids = X - mean(X),
        mean_est = mean(X),
        B = 250,
        block_size = 20,
        seed_addition = i
      )
      se_boot[i] <- results$sd
      covered[i] <- (results$ci_lower < mean_X) & (mean_X < results$ci_upper)
    }

    # Assertions
    expect_lt(
      abs(mean(se_boot) - se_gaussian_expected),
      tol,
      label = paste("Bootstrap SE should approximate theoretical SE for", method)
    )

    expect_true(
      all(se_boot > 0),
      label = paste("All standard errors should be positive for", method)
    )

    coverage <- mean(covered)
    expect_true(
      coverage >= 0.9 && coverage <= 1,
      label = paste("CI coverage should be between 90% and 100% for", method)
    )
  }

  test_bootstrap_method("Bayesian", mean_X, sd_X, n, I, se_gaussian_expected, tol)
  test_bootstrap_method("wild", mean_X, sd_X, n, I, se_gaussian_expected, tol)
  # test_bootstrap_method("naive", mean_X, sd_X, n, I, se_gaussian_expected, tol)
})
