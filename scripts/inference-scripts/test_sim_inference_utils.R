source(here::here("scripts/inference-scripts/0_sim_inference_utils.R"))
needed_libs <- c("tidyverse", "mvtnorm")
devtools::load_all()

test_that("sim_master should work",{
  source(here::here("scripts/inference-scripts/0_sim_inference_utils.R"))
  result <- sim_master(
    iteration = 3,
    N = 600,
    overlap_label = "very_low",
    error_label = "homoskedastic",
    grid_id = 1
  )

  # Check structure
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true("inference_method" %in% names(result))

  # Check required columns
  required_cols <- c("runID", "k", "inference_method", "att_est", "SE",
                     "CI_lower", "CI_upper", "N_T", "ESS_C", "V_E", "V_P",
                     "SATT", "bias")
  expect_true(all(required_cols %in% names(result)))
})

test_that("one_iteration works with pooled inference only", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  # Run with pooled inference only
  result <- one_iteration(
    i = 1,
    k = 4,
    nc = 400,
    nt = 100,
    include_bootstrap = FALSE,
    R = 1,
    verbose = FALSE
  )

  # Check structure
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)  # Only pooled method
  expect_true("inference_method" %in% names(result))
  expect_equal(result$inference_method, "pooled")

  # Check required columns
  required_cols <- c("runID", "k", "inference_method", "att_est", "SE",
                     "CI_lower", "CI_upper", "N_T", "ESS_C", "V_E", "V_P",
                     "SATT", "bias")
  expect_true(all(required_cols %in% names(result)))

  # Check that bootstrap-specific columns are NA for pooled method
  expect_true(is.na(result$boot_mtd))
  expect_true(is.na(result$B))

  # Check that sigma_hat is available for pooled method
  expect_false(is.na(result$sigma_hat))
})



test_that("one_iteration works with different bootstrap methods", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  bootstrap_methods <- c("wild", "Bayesian", "naive")

  for (method in bootstrap_methods) {
    result <- one_iteration(
      i = 1,
      k = 2,
      nc = 50,
      nt = 15,
      include_bootstrap = TRUE,
      boot_mtd = method,
      B = 25,  # Small B for faster testing
      R = 1,
      verbose = FALSE
    )

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 2)

    bootstrap_row <- result[result$inference_method == "bootstrap", ]
    expect_equal(bootstrap_row$boot_mtd, method)

    # Check that we get finite results (or NA if estimation failed)
    expect_true(is.finite(bootstrap_row$att_est) || is.na(bootstrap_row$att_est))
    expect_true(is.finite(bootstrap_row$SE) || is.na(bootstrap_row$SE))
  }
})

test_that("same i returns same result and different i returns different results", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  # Run twice with the same i
  iteration_i <- 2
  set.seed(123 + iteration_i * 10)
  result1 <- one_iteration(
    i = 2,
    k = 4,
    nc = 400,
    nt = 100,
    include_bootstrap = FALSE,
    R = 1,
    verbose = FALSE
  )

  iteration_i <- 2
  set.seed(123 + iteration_i * 10)
  result2 <- one_iteration(
    i = 2,
    k = 4,
    nc = 400,
    nt = 100,
    include_bootstrap = FALSE,
    R = 1,
    verbose = FALSE
  )

  # Run with a different i
  iteration_i <- 3
  set.seed(123 + iteration_i * 10)
  result3 <- one_iteration(
    i = 3,
    k = 4,
    nc = 400,
    nt = 100,
    include_bootstrap = FALSE,
    R = 1,
    verbose = FALSE
  )

  # Check that same i yields identical results
  expect_equal(result1, result2)

  # Check that different i yields different results (in at least one key column)
  expect_false(isTRUE(all.equal(result1$att_est, result3$att_est)))
})


test_that("one_iteration handles different dimensionalities", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  dimensions <- c(1, 2, 4)

  for (k_val in dimensions) {
    result <- one_iteration(
      i = 1,
      k = k_val,
      nc = 80,
      nt = 20,
      include_bootstrap = TRUE,
      B = 25,
      R = 1,
      verbose = FALSE
    )

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 2)
    expect_equal(unique(result$k), k_val)

    # Both methods should have the same sample sizes
    expect_equal(length(unique(result$N_T)), 1)
    expect_equal(length(unique(result$ESS_C)), 1)
  }
})

test_that("one_iteration saves match objects correctly", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  # Create temporary directory
  temp_dir <- tempdir()
  save_path <- file.path(temp_dir, "test_matches")
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

  # Run with save_match = TRUE
  result <- one_iteration(
    i = 999,  # Use unique iteration number
    k = 2,
    nc = 50,
    nt = 15,
    save_match = TRUE,
    path_to_save_match = save_path,
    include_bootstrap = FALSE,
    R = 1,
    verbose = FALSE
  )

  # Check that file was created
  expected_files <- list.files(save_path, pattern = "match_iter_999.*\\.rds$")
  expect_true(length(expected_files) > 0)

  # Check that saved object can be loaded
  if (length(expected_files) > 0) {
    saved_match <- readRDS(file.path(save_path, expected_files[1]))
    expect_true(!is.null(saved_match))
  }

  # Clean up
  unlink(save_path, recursive = TRUE)
})

test_that("one_iteration handles missing data gracefully", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  # Test with very small sample that might cause matching to fail
  result <- suppressWarnings(
    one_iteration(
      i = 1,
      k = 6,  # High dimensionality with small sample
      nc = 10,
      nt = 3,
      include_bootstrap = TRUE,
      B = 10,
      R = 1,
      verbose = FALSE
    )
  )

  expect_true(is.data.frame(result))
  # Should still return results even if estimation fails (with NAs)
  expect_equal(nrow(result), 2)  # Both methods attempted
  expect_true(all(c("pooled", "bootstrap") %in% result$inference_method))
})

test_that("one_iteration reproducibility with seed", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  # Run twice with same seed_addition
  set.seed(123)
  result1 <- one_iteration(
    i = 1,
    k = 2,
    nc = 50,
    nt = 15,
    include_bootstrap = TRUE,
    boot_mtd = "wild",
    B = 25,
    seed_addition = 42,
    R = 1,
    verbose = FALSE
  )

  set.seed(123)
  result2 <- one_iteration(
    i = 1,
    k = 2,
    nc = 50,
    nt = 15,
    include_bootstrap = TRUE,
    boot_mtd = "wild",
    B = 25,
    seed_addition = 42,
    R = 1,
    verbose = FALSE
  )

  # Bootstrap results should be identical with same seed
  bootstrap1 <- result1[result1$inference_method == "bootstrap", ]
  bootstrap2 <- result2[result2$inference_method == "bootstrap", ]

  if (!is.na(bootstrap1$SE) && !is.na(bootstrap2$SE)) {
    expect_equal(bootstrap1$SE, bootstrap2$SE, tolerance = 1e-10)
    expect_equal(bootstrap1$CI_lower, bootstrap2$CI_lower, tolerance = 1e-10)
    expect_equal(bootstrap1$CI_upper, bootstrap2$CI_upper, tolerance = 1e-10)
  }
})

test_that("one_iteration column consistency", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  result <- one_iteration(
    i = 1,
    k = 3,
    nc = 100,
    nt = 25,
    include_bootstrap = TRUE,
    R = 1,
    verbose = FALSE
  )

  # Check that all expected columns are present
  expected_cols <- c(
    "runID", "k", "inference_method", "boot_mtd", "B",
    "att_est", "SE", "CI_lower", "CI_upper",
    "N_T", "ESS_C", "V_E", "V_P", "sigma_hat",
    "SATT", "bias"
  )

  expect_true(all(expected_cols %in% names(result)))

  # Check column types
  expect_true(is.numeric(result$runID))
  expect_true(is.numeric(result$k))
  expect_true(is.character(result$inference_method))
  expect_true(is.numeric(result$att_est))
  expect_true(is.numeric(result$SE))
  expect_true(is.numeric(result$SATT))
})

# Integration test that mimics the actual usage
test_that("Integration test: Full workflow", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  # Test multiple iterations (small scale)
  results_list <- list()

  for (i in 1:2) {
    results_list[[i]] <- one_iteration(
      i = i,
      k = 2,
      nc = 80,
      nt = 20,
      include_bootstrap = TRUE,
      boot_mtd = "wild",
      B = 30,
      R = 2,
      verbose = FALSE
    )
  }

  # Combine results
  all_results <- bind_rows(results_list)

  expect_equal(nrow(all_results), 4)  # 2 iterations × 2 methods
  expect_equal(length(unique(all_results$runID)), 2)
  expect_true(all(c("pooled", "bootstrap") %in% all_results$inference_method))

  # Check that we can separate and analyze results by method
  pooled_results <- all_results %>% filter(inference_method == "pooled")
  bootstrap_results <- all_results %>% filter(inference_method == "bootstrap")

  expect_equal(nrow(pooled_results), 2)
  expect_equal(nrow(bootstrap_results), 2)
})
