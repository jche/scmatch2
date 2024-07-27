


# functions for estimating effects
# The following two functions should be combined
#   once we get to know the input of get_att_ests
get_est_att_from_wt <- function(df,
                                input_wt){
  df_est_att <- df %>%
    cbind(wt=input_wt) %>%
    group_by(Z) %>%
    summarise(Y_wtd = weighted.mean(Y,wt))

  est_att <- diff(df_est_att$Y_wtd)
  return(est_att)
}


#' Estimate the ATT from a matched dataframe
#'
#' Given a matched dataset, calculate the estimated ATT
#'
#' @param matched_df A matched dataset
#'
#' @export
get_att_ests <- function(matched_df) {
  if ( is.csm_matches(matched_df) ) {
    matched_df <- matched_df$result
  }
  matched_df %>%
    group_by(Z) %>%
    summarize(mn = sum(Y*weights) / sum(weights)) %>%
    summarize(est = last(mn) - first(mn)) %>%
    pull(est)
}

# get_att_ests_from_wt<-
#   function(df, input_wt)
