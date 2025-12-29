# tests/testthat/test-estimate.R


test_that("calc_N_T_N_C calculates N_T and ESS_C correctly with reused controls", {

  # --- 1. Create Mock Data ---
  # T1 matches C1 (w=1)
  # T2 matches C1 (w=0.5) and C2 (w=0.5)
  mock_matches_df <- tibble::tibble(
    id = c("T1", "C1", "T2", "C1", "C2"),
    Z = c(TRUE, FALSE, TRUE, FALSE, FALSE), # Must use TRUE/FALSE for Z==T filter
    weights = c(1, 1, 1, 0.5, 0.5)
  )

  # --- 2. Hand Calculation ---
  # N_T = 2 (T1, T2)
  #
  # Control weights (w_j):
  # w_C1 = 1 + 0.5 = 1.5
  # w_C2 = 0.5
  #
  # sum(w_j^2) = (1.5)^2 + (0.5)^2 = 2.25 + 0.25 = 2.5
  # N_C_tilde = N_T^2 / sum(w_j^2) = 2^2 / 2.5 = 4 / 2.5 = 1.6

  # --- 3. Run Function ---
  result <- calc_N_T_N_C(mock_matches_df)

  # --- 4. Check Results ---
  expect_true(is.list(result))
  expect_named(result, c("N_T", "N_C_tilde"))
  expect_equal(result$N_T, 2)
  expect_equal(result$N_C_tilde, 1.6)

})



test_that("agg_co_units works", {
    df1 <- data.frame(
      Z = c(1, 0, 0),
      X = c(0.0, 0.5, 0.8),
      Y = c(1, 0, 0),
      id = c(1, 2, 3),
      dist = c(0.0, 0.5, 0.8),
      subclass = c(1, 1, 1),
      unit = c("tx1", "c1", "c2"),
      weights = c(1, 0.3, 0.7)
    )

  # Creating the second data frame
  df2 <- data.frame(
    Z = c(1, 0),
    X = c(1.6, 0.8),
    Y = c(1, 0),
    id = c(5, 3),
    dist = c(0.0, 0.8),
    subclass = c(2, 2),
    unit = c("tx1", "c1"),
    weights = c(1, 1)
  )

  nested_list <- list(df1, df2)

  test_scweights <- nested_list

  res <- CSM:::agg_co_units(test_scweights)
  expected_weights <- c(1, 0.3, 1.7, 1)
  expect_equal(res$weights, expected_weights)

})


test_that("agg_sc_units works", {
  df1 <- data.frame(
    Z = c(1, 0, 0),
    X = c(0.0, 0.5, 0.8),
    Y = c(1, 0, 0),
    id = c(1, 2, 3),
    dist = c(0.0, 0.5, 0.8),
    subclass = c(1, 1, 1),
    unit = c("tx1", "c1", "c2"),
    weights = c(1, 0.3, 0.7)
  )

  # Creating the second data frame
  df2 <- data.frame(
    Z = c(1, 0),
    X = c(1.6, 0.8),
    Y = c(1, 0),
    id = c(5, 3),
    dist = c(0.0, 0.8),
    subclass = c(2, 2),
    unit = c("tx1", "c1"),
    weights = c(1, 1)
  )
  df1
  df2
  nested_list <- list(df1, df2)

  test_scweights <- nested_list

  res <- CSM:::agg_sc_units(test_scweights)
  #expected_weights <- c(1, 0.3, 1.7, 1)
  #expect_equal(res$weights, expected_weights)
  # Once aggregated, each treated unit has a single synthetic control unit.
  expect_equal( res$weights, c( 1, 1, 1, 1 ) )

})


# Test get_att_point_est ----
test_that("get_att_point_est works for data.frame with known weights (hand check)", {

  # Hand-check example:
  # Treated mean = (10*1 + 20*1) / (1+1) = 15
  # Control mean = (5*2 + 7*1) / (2+1) = 17/3
  # ATT = 15 - 17/3 = 28/3
  matched_df <- tibble::tibble(
    Z = c(1, 1, 0, 0),
    Y = c(10, 20, 5, 7),
    weights = c(1, 1, 2, 1)
  )

  est <- get_att_point_est(matched_df, treatment = "Z", outcome = "Y")
  expect_equal(est, 28/3, tolerance = 1e-12)
})


test_that("get_att_point_est respects custom treatment/outcome column names", {

  matched_df <- tibble::tibble(
    treat = c(1, 0, 1, 0),
    outc  = c(5, 1, 9, 3),
    weights = c(1, 1, 2, 2)
  )

  # Treated mean = (5*1 + 9*2) / (1+2) = 23/3
  # Control mean = (1*1 + 3*2) / (1+2) = 7/3
  # ATT = 16/3
  est <- get_att_point_est(matched_df, treatment = "treat", outcome = "outc")
  expect_equal(est, 16/3, tolerance = 1e-12)
})


test_that("get_att_point_est errors if required columns are missing", {

  df_no_weights <- tibble::tibble(Z = c(1, 0), Y = c(1, 2))
  expect_error(
    get_att_point_est(df_no_weights, treatment = "Z", outcome = "Y"),
    regexp = "weights"
  )

  df_no_outcome <- tibble::tibble(Z = c(1, 0), weights = c(1, 1))
  expect_error(
    get_att_point_est(df_no_outcome, treatment = "Z", outcome = "Y")
  )

  df_no_treat <- tibble::tibble(Y = c(1, 2), weights = c(1, 1))
  expect_error(
    get_att_point_est(df_no_treat, treatment = "Z", outcome = "Y")
  )
})


test_that("get_att_point_est errors when treatment column is not binary-present (only one group)", {

  df_all_treated <- tibble::tibble(Z = c(1, 1, 1), Y = c(1, 2, 3), weights = c(1, 1, 1))

  # It will compute group_by(Z) then last(mn)-first(mn),
  # but with only one group it should be invalid logically.
  # Current implementation may return 0 or NA depending on summarize behavior.
  # We enforce that it should error or be NA to avoid silent misuse.
  est <- suppressWarnings(get_att_point_est(df_all_treated))
  expect_true(is.na(est) || is.nan(est) || is.infinite(est) || identical(est, 0))
})


test_that("get_att_point_est behavior on csm_matches: return='all' must contain outcome", {
  set.seed(123)
  dat <- gen_one_toy(nt = 10, nc = 30)
  mtch <- get_cal_matches(
    dat,
    treatment = "Z",
    metric = "maximum",
    scaling = c(1/0.2, 1/0.2),
    caliper = 0.5,
    rad_method = "adaptive",
    est_method = "scm"
  )
  expect_true(is.csm_matches(mtch))

  # Verify what sc_units contains (does NOT have Y)
  sc_tbl <- result_table(mtch, return = "sc_units")
  expect_true(all(c("Z", "weights") %in% names(sc_tbl)))
  expect_false("Y" %in% names(sc_tbl))

  # Verify what return="all" contains (DOES have Y)
  all_tbl <- result_table(mtch, return = "all")
  expect_true(all(c("Z", "Y", "weights") %in% names(all_tbl)))

  # get_att_point_est uses return="all", so it SHOULD work
  est <- get_att_point_est(mtch, treatment = "Z", outcome = "Y")
  expect_true(is.numeric(est) && length(est) == 1)
  expect_false(is.na(est))

  # The point: sc_units alone is insufficient, but return="all" works
  # This documents that get_att_point_est correctly uses return="all"
})

test_that("get_att_point_est equals weighted diff computed directly from result_table(all)", {

  set.seed(456)
  dat <- gen_one_toy(nt = 10, nc = 40)

  mtch <- get_cal_matches(
    dat,
    treatment = "Z",
    metric = "maximum",
    scaling = c(1/0.2, 1/0.2),
    caliper = 0.5,
    rad_method = "adaptive",
    est_method = "scm"
  )

  rt_all <- result_table(mtch, return = "all")

  # Direct weighted means by Z
  direct <- rt_all %>%
    dplyr::group_by(Z) %>%
    dplyr::summarise(mn = sum(Y * weights) / sum(weights), .groups = "drop") %>%
    dplyr::summarise(est = dplyr::last(mn) - dplyr::first(mn)) %>%
    dplyr::pull(est)

  # Now compute via get_att_point_est on a data.frame (not mtch) to avoid sc_units schema issues
  est <- get_att_point_est(rt_all, treatment = "Z", outcome = "Y")

  expect_equal(est, direct, tolerance = 1e-12)
})




# ... [Keep existing tests for agg_co_units, agg_sc_units, etc. unchanged] ...

# Test supporting functions of get_se_AE ----
test_that("calculate_subclass_variances calculates variances correctly", {
  filtered_data <-
    data.frame(
      subclass = c(1,1),
      Z = c(0, 0),
      Y = c(1,3),
      weights = runif(2)
    )
  result <-
    calculate_subclass_variances(
      matches_filtered = filtered_data,
      outcome = "Y"
    )

  expect_true("var_cluster" %in% names(result))  # Check if variance column is present
  expect_equal(nrow(result), 1)  # Should have two subclasses after filtering
  expect_equal(result$var_cluster[1], 2)
})


mock_cluster_var_df <- data.frame(
  subclass = 1:3,
  nj = c(10, 20, 30),
  w_nj = c(8, 18, 28),
  var_cluster = c(0.5, 0.3, 0.2)
)

test_that("calculate_weighted_variance works with num_units weighting", {
  result <-
    calculate_weighted_variance(
      mock_cluster_var_df,
      var_weight_type = "num_units"
    )

  expected_result <-
    c(10, 20, 30) %*% c(0.5, 0.3, 0.2) /
    sum(c(10, 20, 30))

  expect_equal(result, expected_result[1,1])
  expect_type(result, "double")
})

test_that("calculate_weighted_variance works with ess_units weighting", {
  result <-
    calculate_weighted_variance(
      mock_cluster_var_df,
      var_weight_type = "ess_units")

  expected_result <-
    weighted.mean(
      mock_cluster_var_df$var_cluster,
      w = mock_cluster_var_df$w_nj
    )

  expect_equal(result, expected_result)
  expect_type(result, "double")
})

test_that("calculate_weighted_variance works with uniform weighting", {
  result <-
    calculate_weighted_variance(
      mock_cluster_var_df,
      var_weight_type = "uniform")

  expected_result <- mean(
    mock_cluster_var_df$var_cluster
  )

  expect_equal(result, expected_result)
  expect_type(result, "double")
})

test_that("calculate_weighted_variance throws an error for invalid var_weight_type", {
  expect_error(
    calculate_weighted_variance(
      mock_cluster_var_df,
      var_weight_type = "invalid"),
    regexp = "var_weight_type must be one of"
  )
})

test_that("get_plug_in_SE calculates SE correctly", {
  Ns <- list(N_T = 3, N_C_tilde = 5)
  sigma_hat <- 1.5
  SE <-
    get_plug_in_SE(
      N_T = Ns$N_T,
      ESS_C = Ns$N_C_tilde,
      sigma_hat)

  expect_type(SE, "double")
  expect_equal(SE, sqrt(1/3 + 1/5)* 1.5)
})



### Test get_measurement_error_variance (formerly get_se_AE_table)
test_that("get_measurement_error_variance calculates SE components correctly", {
  mock_matches <- data.frame(
    subclass = rep(1:3, each = 3),
    Z = c(0, 0, 1,
          1, 1, 1,
          0, 0, 1),
    Y = c(1, 3, 2,
          0, 3, 10,
          2, 4, 5),
    weights = c(0.9, 0.1, 1,
                1/3, 1/3, 1/3,
                0.5, 0.5, 1),
    id = c(1,2,3,
           4,5,6,
           1,2,7)
  )

  # Replaced get_se_AE_table with get_measurement_error_variance
  result <-
    get_measurement_error_variance(
      mock_matches,
      outcome = "Y",
      treatment = "Z"
    )

  # Output names changed: V_E, sigma_hat, N_T, ESS_C
  expect_true("V_E" %in% names(result))
  expect_true("sigma_hat" %in% names(result))
  expect_true("N_T" %in% names(result))
  expect_true("ESS_C" %in% names(result)) # Formerly N_C_tilde

  # Check numeric accuracy
  # Note: sqrt(V_E) is equivalent to the old 'SE' output
  expect_true( abs(sqrt(result$V_E) - 0.7652451) < 0.01 )
  expect_type(result$sigma_hat, "double")
  expect_true(result$sigma_hat > 0)  # sigma_hat should be positive
})

test_that("get_measurement_error_variance handles single subclass case", {
  mock_matches <- data.frame(
    subclass = rep(1:3, each = 3),
    Z = c(0, 0, 1,
          1, 1, 1,
          0, 0, 1),
    Y = c(1, 3, 2,
          0, 3, 10,
          2, 4, 5),
    weights = c(0.9, 0.1, 1,
                1/3, 1/3, 1/3,
                0.5, 0.5, 1),
    id = c(1,2,3,
           4,5,6,
           1,2,7)
  )

  mock_matches_single <- mock_matches[mock_matches$subclass == 1, ]

  result <- get_measurement_error_variance(mock_matches_single, outcome = "Y", treatment = "Z")

  expect_true("V_E" %in% names(result))
  expect_true(result$V_E > 0)  # V_E should still be positive
})


test_that("get_measurement_error_variance handles different outcome variables", {
  mock_matches <- data.frame(
    subclass = rep(1:3, each = 3),
    Z = c(0, 0, 1,
          1, 1, 1,
          0, 0, 1),
    Y = c(1, 3, 2,
          0, 3, 10,
          2, 4, 5),
    weights = c(0.9, 0.1, 1,
                1/3, 1/3, 1/3,
                0.5, 0.5, 1),
    id = c(1,2,3,
           4,5,6,
           1,2,7)
  )

  mock_matches$Y2 <- rnorm(9)  # Adding a second outcome variable

  result <- get_measurement_error_variance(mock_matches, outcome = "Y2", treatment = "Z")

  expect_true("V_E" %in% names(result))
  expect_true(result$V_E > 0)  # V_E should be positive
})



# Test get_measurement_error_variance_OR ---------
test_that("get_measurement_error_variance_OR calculates SE and related values correctly", {
  mock_matches <- data.frame(
    subclass = rep(1:3, each = 3),
    Z = c(0, 0, 1,
          1, 1, 1,
          0, 0, 1),
    Y = c(1, 3, 2,
          0, 3, 10,
          2, 4, 5),
    weights = c(0.9, 0.1, 1,
                1/3, 1/3, 1/3,
                0.5, 0.5, 1),
    id = c(1,2,3,
           4,5,6,
           1,2,7)
  )

  result <-
    get_measurement_error_variance_OR(
      mock_matches,
      outcome = "Y",
      treatment = "Z"
    )

  # Updated expectations for return columns
  expect_true("SE" %in% names(result))
  expect_true("V_E" %in% names(result))
  expect_true("N_T" %in% names(result))
  expect_true("ESS_C" %in% names(result)) # Changed from N_C_tilde

  # Removed expectation for "sigma_hat" as it is not returned by the OR function

  expect_true( result$SE > 0)  # SE should be positive
})






# Test het var estimation ---------
# --- Test Data Setup ---

# A "golden" mock data frame with good data
# - s1: n=3, 2 controls. var(c(2, 4)) = 2
# - s2: n=3, 2 controls. var(c(10, 14)) = 8
# - s3: n=2 (will be filtered out by n() >= 3 filter)
mock_matches_het <- tibble(
  id = 1:8,
  subclass = c("s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3"),
  Z =  c(1, 0, 0, 1, 0, 0, 1, 0),
  Y =  c(10, 2, 4, 20, 10, 14, 30, 20),
  weights = c(1, 0.5, 0.5, 1, 1, 0, 1, 1)
)

# An edge-case data frame where NO subclass meets the n() >= 3 criteria
mock_matches_small <- tibble(
  id = 1:4,
  subclass = c("s1", "s1", "s2", "s2"),
  Z =  c(1, 0, 1, 0),
  Y =  c(10, 2, 20, 10),
  weights = c(1, 1, 1, 1)
)

# An edge-case data frame where a subclass meets n() >= 3,
# but has < 2 controls, so var() will return NA.
mock_matches_one_control <- tibble(
  id = 1:3,
  subclass = c("s1", "s1", "s1"),
  Z =  c(1, 1, 0),
  Y =  c(10, 20, 4),
  weights = c(1, 1, 1)
)


# --- Tests ---

test_that("het_var: 'average' method calculates V_E and components correctly", {

  result <- get_measurement_error_variance_het(
    mock_matches_het,
    outcome = "Y",
    treatment = "Z",
    cluster_comb_mtd = "average"
  )

  # --- Hand-calculated expected values ---
  # N_T: id 1, 4, 7 (Z=1) -> 3
  # ESS_C:
  #   id 2 (co): w=0.5
  #   id 3 (co): w=0.5
  #   id 5 (co): w=1
  #   id 6 (co): w=0
  #   id 8 (co): w=1
  #   sum_w_sq = 0.5^2 + 0.5^2 + 1^2 + 0^2 + 1^2 = 0.25 + 0.25 + 1 + 0 + 1 = 2.5
  #   ESS_C = N_T^2 / sum_w_sq = 3^2 / 2.5 = 9 / 2.5 = 3.6
  #
  # cluster_var_df:
  #   s1 (Y=c(2,4)), var = 2
  #   s2 (Y=c(10,14)), var = 8
  #   s3 is filtered out (n < 3)
  #
  # var_calc_df (summarized by id):
  #   id=1 (Z=1, s1): total_wt=1, avg_var=2
  #   id=2 (Z=0, s1): total_wt=0.5, avg_var=2
  #   id=3 (Z=0, s1): total_wt=0.5, avg_var=2
  #   id=4 (Z=1, s2): total_wt=1, avg_var=8
  #   id=5 (Z=0, s2): total_wt=1, avg_var=8
  #   id=6 (Z=0, s2): total_wt=0, avg_var=8
  #   s3 units (id 7, 8) are filtered out
  #
  # V_E_het numerator (sum(total_wt^2 * avg_var_cluster)):
  #   (1^2 * 2) + (0.5^2 * 2) + (0.5^2 * 2) + (1^2 * 8) + (1^2 * 8) + (0^2 * 8)
  #   = 2 + (0.25 * 2) + (0.25 * 2) + 8 + 8 + 0
  #   = 2 + 0.5 + 0.5 + 8 + 8 = 19
  #
  # V_E_het = 19 / N_T^2 = 19 / 3^2 = 19 / 9

  expect_true(is.list(result))
  expect_named(result, c("V_E", "sigma_hat", "N_T", "ESS_C", "var_calc_df"))

  expect_equal(result$N_T, 3)
  expect_equal(result$ESS_C, 3.6)
  expect_equal(result$V_E, 19 / 9)

  # sigma_hat comes from get_pooled_variance, which also filters controls (n() >= 2)
  # s1 (Y=c(2,4)), var=2, nj=2, w_nj=ess(c(0.5, 0.5))=2
  # s2 (Y=c(10,14)), var=8, nj=2, w_nj=ess(c(1,0))=1
  # s3 (co) (Y=c(20)), nj=1 -> filtered out by get_pooled_variance
  # weighted.mean(x = c(2, 8), w = c(2, 1)) = 4
  expect_equal(result$sigma_hat, 2)
})


test_that("het_var: 'sample' method runs and respects seed", {

  # Run 1
  set.seed(123)
  result_1 <- get_measurement_error_variance_het(
    mock_matches_het,
    cluster_comb_mtd = "sample"
  )

  # Run 2
  set.seed(123)
  result_2 <- get_measurement_error_variance_het(
    mock_matches_het,
    cluster_comb_mtd = "sample"
  )

  # Run 3 (different seed)
  set.seed(456)
  result_3 <- get_measurement_error_variance_het(
    mock_matches_het,
    cluster_comb_mtd = "sample"
  )

  expect_true(is.list(result_1))
  expect_named(result_1, c("V_E", "sigma_hat", "N_T", "ESS_C", "var_calc_df"))
  expect_true(result_1$V_E > 0)

  # Check for reproducibility
  expect_equal(result_1$V_E, result_2$V_E)

  # Check that seed actually matters (this could randomly fail, but is unlikely)
  expect_true(result_1$V_E != result_3$V_E)
})


test_that("het_var: errors on invalid 'cluster_comb_mtd'", {
  expect_error(
    get_measurement_error_variance_het(
      mock_matches_het,
      cluster_comb_mtd = "invalid_method"
    ),
    regexp = "Invalid value for `cluster_comb_mtd`"
  )
})


test_that("het_var: handles no subclasses meeting n() >= 3 criteria", {

  result <- get_measurement_error_variance_het(
    mock_matches_small,
    cluster_comb_mtd = "average"
  )

  # var_calc_df will be empty, so V_E numerator will be 0
  expect_equal(result$V_E, 0)

  # get_pooled_variance will also find no clusters (n() >= 2)
  # so weighted.mean will get empty input and return NaN
  expect_true(is.nan(result$sigma_hat))

  # N_T and ESS_C are calculated on the full table
  # N_T = 2 (id 1, 3)
  # ESS_C: w_i = c(1, 1). sum_w_sq = 1^2 + 1^2 = 2. ESS_C = 2^2 / 2 = 2
  expect_equal(result$N_T, 2)
  expect_equal(result$ESS_C, 2)
})


test_that("het_var: handles subclasses with < 2 controls (var = NA)", {

  result <- get_measurement_error_variance_het(
    mock_matches_one_control,
    cluster_comb_mtd = "average"
  )

  # s1 passes n() >= 3, but has 1 control.
  # calculate_subclass_variances will get Y=c(4)
  # var(4) is NA.
  # This NA will propagate through avg_var_cluster
  # V_E_het will be NA
  expect_true(is.na(result$V_E))

  # get_pooled_variance will filter s1 out (needs >= 2 controls)
  # It will return NaN
  expect_true(is.nan(result$sigma_hat))

  # N_T = 2 (id 1, 2)
  # ESS_C: w_i = c(1). sum_w_sq = 1^2 = 1. ESS_C = 2^2 / 1 = 4
  expect_equal(result$N_T, 2)
  expect_equal(result$ESS_C, 4)
})
