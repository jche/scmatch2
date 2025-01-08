
# Test create_love_plot_df
# Conclusion:
#   - return = "sc_units" is a must
#   - return = "agg_co_units" cannot make it work

test_that("diagnostic plots", {

  test_df <-
    tibble(Z=c(1,0,0,0,1),
           X=c(0,-0.5,0.8,3,1.6) ) %>%
    mutate(
      X2 = rnorm(n()),
      X3 = rnorm(n()) + Z )

  res <- get_cal_matches(df = test_df,
                         covs = "X",
                         treatment = "Z",
                         metric = "maximum",
                         caliper = 1,
                         rad_method = "adaptive",
                         est_method = "scm",
                         scaling = 1)

  feasible <- feasible_units(res)
  n_feasible <- nrow(feasible)

  # It seems telling me that the NA is the aggregated control so it's one unit
  # attr(res, "adacalipers")
  covs = c( "X", "X2","X3" )

  love_plot_df <-
    CSM:::create_love_plot_df( res = res,
                               covs = covs)
  love_plot_df

  expect_true( is.data.frame(love_plot_df) )

} )
