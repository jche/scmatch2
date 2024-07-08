test_that("get_SL_pred returns correct output structure", {
  mock_SL_fit <- create_mock_SL_fit()
  mock_df_test <- create_mock_df_test()
  X_names <- c("X1","X2")
  result <- get_SL_pred(mock_SL_fit,
                        mock_df_test,
                        X_names)

  expect_true(is.matrix(result) == T)
})



test_that("get_SL_fit returns a correct output data type", {
  result = create_mock_SL_fit()
  expect_true("SuperLearner" %in% class(result))
})
