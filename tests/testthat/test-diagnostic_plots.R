
# Test create_love_plot_df
# Conclusion:
#   - return = "sc_units" is a must
#   - return = "agg_co_units" cannot make it work
test_df <-
  data.frame(Z=c(1,0,0,0,1),
             X=c(0,-0.5,0.8,3,1.6))
attr(res, "feasible_units")
res <- get_cal_matches(df = test_df,
                       covs = "X",
                       treatment = "Z",
                       metric = "maximum",
                       caliper = 1,
                       rad_method = "adaptive",
                       est_method = "scm",
                       return = "sc_units",
                       # return = "agg_co_units",
                       dist_scaling = 1)
feasible_subclasses <-
  attr(res, "feasible_subclasses")
n_feasible <- length(feasible_subclasses)

# It seems telling me that the NA is the aggregated control so it's one unit
# attr(res, "adacalipers")
covs = "X"

love_plot_df <-
  create_love_plot_df(res = res,
                      covs = "X")
