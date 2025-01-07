
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

######
# Test supporting functions of get_se_AE
######

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

######
# Test get_se_AE
######
test_that("get_se_AE_table calculates SE and related values correctly", {
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
    get_se_AE_table(
      mock_matches,
      outcome = "Y",
      treatment = "Z"
    )

  expect_true("SE" %in% names(result))
  expect_true("sigma_hat" %in% names(result))
  expect_true("N_T" %in% names(result))
  expect_true("N_C_tilde" %in% names(result))

  expect_true( abs(result$SE - 0.7652451) < 0.01 )  # SE should be positive
  expect_type(result$sigma_hat, "double")
  expect_true(result$sigma_hat > 0)  # sigma_hat should be positive
})

test_that("get_se_AE_table handles single subclass case", {
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

  result <- get_se_AE_table(mock_matches_single, outcome = "Y", treatment = "Z")

  expect_true("SE" %in% names(result))
  expect_true(result$SE > 0)  # SE should still be positive
})


test_that("get_se_AE_table handles different outcome variables", {
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

  result <- get_se_AE_table(mock_matches, outcome = "Y2", treatment = "Z")

  expect_true("SE" %in% names(result))
  expect_true(result$SE > 0)  # SE should be positive
})
