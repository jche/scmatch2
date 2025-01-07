
test_that("Bootstrap standard errors match theoretical expectations for Gaussian data", {
  set.seed(42)
  mean_X <- 0
  sd_X <- 2
  n <- 500
  I <- 100
  se_gaussian_expected <- sd_X / sqrt(n)
  tol <- 0.01

  test_bootstrap_method <- function(method, mean_X, sd_X, n, I, se_gaussian_expected, tol) {
    # Create bootstrap CI calculator for the given method
    bootstrap_ci <- make_bootstrap_ci(method)

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
  # test_bootstrap_method("naive-resid", mean_X, sd_X, n, I, se_gaussian_expected, tol)
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
  # test_bootstrap_method("naive-resid", mean_X, sd_X, n, I, se_gaussian_expected, tol)
})
