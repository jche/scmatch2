test_that("gen_toy_covar returns correct output",{
  n = 100
  X1_ctrs <- c(0.25, 0.75)
  X2_ctrs <- c(0.25, 0.75)
  SD <- 0.1
  res <- gen_toy_covar(n, X1_ctrs, X2_ctrs, SD)

  expect_true(is_tibble(res))

  expect_true(all(c("X1", "X2") %in% names(res) ))

  expect_equal(nrow(res),n)

  expect_equal(mean(res$X1), mean(X1_ctrs),
               tolerance = 0.1)

  expect_equal(mean(res$X2), mean(X2_ctrs),
               tolerance = 0.1)
  expect_equal(sd(res$X1), SD, tolerance = 0.5)
  expect_equal(sd(res$X2), SD, tolerance = 0.5)
})
