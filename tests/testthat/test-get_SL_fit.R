library(CSM)
library(tidyverse)

test_that("get_SL_fit returns a correct output data type", {
  result = create_mock_SL_fit()
  expect_true("SuperLearner" %in% class(result))
})
