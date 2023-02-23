
# functions for estimating treatment effects


get_att_ests <- function(matched_df) {
  matched_df %>%
    group_by(Z) %>%
    summarize(mn = sum(Y*weights) / sum(weights)) %>%
    summarize(est = last(mn) - first(mn)) %>%
    pull(est)
}

