source("./scripts/inference-scripts/0_sim_inference_utils.R")
needed_libs <- c("tidyverse", "mvtnorm")

test_that("one_iteration works with pooled inference only", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  # Run with pooled inference only
  result <- one_iteration(
    i = 1,
    k = 2,
    nc = 100,
    nt = 25,
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

test_that("one_iteration works with both pooled and bootstrap inference", {
  skip_if_not(load_libs(needed_libs), "Required libraries not available")

  # Run with both methods
  result <- one_iteration(
    i = 1,
    k = 2,
    nc = 100,
    nt = 25,
    include_bootstrap = TRUE,
    boot_mtd = "wild",
    B = 50,  # Small B for faster testing
    R = 1,
    verbose = FALSE
  )

  # Check structure
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)  # Both pooled and bootstrap
  expect_true(all(c("pooled", "bootstrap") %in% result$inference_method))

  # Check that both rows have the same basic info
  expect_equal(length(unique(result$runID)), 1)
  expect_equal(length(unique(result$k)), 1)
  expect_equal(length(unique(result$SATT)), 1)
  expect_equal(length(unique(result$bias)), 1)

  # Check bootstrap-specific columns
  bootstrap_row <- result[result$inference_method == "bootstrap", ]
  expect_equal(bootstrap_row$boot_mtd, "wild")
  expect_equal(bootstrap_row$B, 50)

  pooled_row <- result[result$inference_method == "pooled", ]
  expect_true(is.na(pooled_row$boot_mtd))
  expect_true(is.na(pooled_row$B))

  # ATT estimates should be the same (point estimate doesn't change)
  if (!is.na(bootstrap_row$att_est) && !is.na(pooled_row$att_est)) {
    expect_equal(bootstrap_row$att_est, pooled_row$att_est, tolerance = 1e-10)
  }
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


# if ( FALSE ) {
#   # --- Test run of one_iteration ---
#   message("--- Testing one_iteration ---")
#   needed_libs_iter <- c("tidyverse", "mvtnorm") # Add CSM if needed globally
#   if (load_libs(needed_libs_iter)) {
#     # Assuming necessary functions (gen_one_toy, etc.) are loaded
#     # You might need to explicitly load CSM or its functions if not attached
#     # e.g. library(CSM) or source("path/to/CSM_functions.R")
#
#     source("./scripts/inference-scripts/0_sim_inference_utils.R")
#     message("Running one_iteration with k=4...")
#     iteration_result_k4 <- one_iteration( i = 1, k = 4, nc = 350, nt = 30, R = 1, verbose=TRUE )
#     print(iteration_result_k4)
#
#     iteration_result_k4 <- one_iteration( i = 1, k = 4, nc = 500, scaling = 8, toy_ctr_dist = 0.5, R = 1, verbose=TRUE )
#     print(iteration_result_k4)
#
#     message("\nRunning one_iteration with k=2...")
#     source("./scripts/inference-scripts/0_sim_inference_utils.R")
#     iteration_result_k2 <- one_iteration( i = 1, k = 2, nc = 100, nt = 25, R = 1, verbose=TRUE )
#     print(iteration_result_k2)
#   } else {
#     message("Skipping one_iteration test due to missing packages.")
#   }
# }
#
# if ( FALSE ) {
#   # --- Test the save_match feature ---
#   message("--- Testing save_match ---")
#   needed_libs_iter <- c("tidyverse", "mvtnorm") # Add CSM if needed globally
#   if (load_libs(needed_libs_iter)) {
#     source("./scripts/inference-scripts/0_sim_inference_utils.R")
#     iteration_result_k4 <- one_iteration( i = 1, k = 4, nc = 350, nt = 30, R = 1,
#                                           save_match = TRUE,
#                                           path_to_save_match = "./data/outputs/4d_bias_inference_Apr2025/test",
#                                           verbose=TRUE )
#   } else {
#     message("Skipping one_iteration test due to missing packages.")
#   }
# }
#
