# Driver functions for core CSM method



get_att_csm <- function(df,
                        covs = starts_with("X"),
                        treatment = "Z",
                        metric = c("maximum", "euclidean", "manhattan"),
                        caliper = 1,
                        rad_method = c("adaptive", "fixed", "1nn"),
                        est_method = c("scm", "scm_extrap", "average"),
                        return = c("sc_units", "agg_co_units", "all"),
                        dist_scaling = df %>%
                          summarize(across(starts_with("X"),
                                           function(x) {
                                             if (is.numeric(x)) 1/sd(x)
                                             else 1000
                                           })),
                        ...) {
  match_weighted_df <- df %>%
    get_cal_matches(
      covs = covs,
      treatment = treatment,
      metric = metric,
      caliper = caliper,
      rad_method = rad_method,
      est_method = est_method,
      return = return,
      dist_scaling = dist_scaling,
      ...
    )
  est <- get_att_ests(match_weighted_df)
  return(est)
}
