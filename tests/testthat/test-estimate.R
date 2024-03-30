
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

  res <- agg_co_units(test_scweights)
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

  nested_list <- list(df1, df2)

  test_scweights <- nested_list

  res <- agg_sc_units(test_scweights)
  # expected_weights <- c(1, 0.3, 1.7, 1)
  #
  # expect_equal(res$weights, expected_weights)
})



