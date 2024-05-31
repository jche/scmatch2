test_that("fit_CSM works well",{
  test_df <-
    data.frame(Z=c(1,0,0,0,1),
               X=c(0,0.5,0.8,3,1.6),
               Y = c(1,0,0,0,1))
  res <- CSM:::fit_CSM(df = test_df,
                       covs = "X",
                       treatment = "Z",
                       metric = "maximum",
                       caliper = 1,
                       rad_method = "adaptive",
                       est_method = "scm",
                       return = "agg_co_units",
                       dist_scaling = 1)
  expect_equal(res, 1)
})
