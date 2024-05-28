
# functions for estimating effects
# The following two functions should be combined
#   once we get to know the input of get_att_ests
get_est_att_from_wt <-
  function(df,
           input_wt){
    df_est_att <- df %>%
      cbind(wt=input_wt) %>%
      group_by(Z) %>%
      summarise(Y_wtd = weighted.mean(Y,wt))

    est_att <- diff(df_est_att$Y_wtd)
    return(est_att)
  }



get_att_ests <- function(matched_df) {
  matched_df %>%
    group_by(Z) %>%
    summarize(mn = sum(Y*weights) / sum(weights)) %>%
    summarize(est = last(mn) - first(mn)) %>%
    pull(est)
}

# get_att_ests_from_wt<-
#   function(df, input_wt)

# aggregation functions ---------------------------------------------------


#' aggregate scweights to tx/sc units
#'
#' Aggregate by cluster defined by treated unit, calculating the
#' weighted average of the control units in the cluster for all
#' variables that start with "X" or "Y" in the dataframes.
#'
#' @param scweights list of sc weights (tibbles with tx unit in first
#'   row).
#'
#' @return df with tx and corresponding sc units
agg_sc_units <- function(scweights) {

  if (!is.data.frame(scweights)) {
    scweights <- scweights %>%
      map_dfr(~mutate(., subclass=id[1]))
  }

  rs <- scweights %>%
    group_by(subclass, Z) %>%
    summarize(across(starts_with("X"),
                     ~sum(.x * weights)),
              across(starts_with("Y"),
                     ~sum(.x * weights)),
              .groups="drop_last") %>%
    mutate(id = c(NA, subclass[1]), .before="subclass") %>%
    mutate(weights = 1) %>%
    ungroup()

  # Make IDs for the synthetic controls
  rs$id = as.character(rs$id)
  rs$id[ is.na(rs$id) ] <- paste0( rs$subclass[ is.na(rs$id) ], "_syn" )

  rs
}


#' aggregate sc weights to original tx and (weighted) co units
#'
#' @param scweights list of sc weights (tibbles with tx unit in first row)
#'
#' @return df with original tx and sc units, with weights
agg_co_units <- function(scweights) {
  if (!is.data.frame(scweights)) {
    scweights <- scweights %>%
      bind_rows()
  }

  scweights %>%
    group_by(id) %>%
    summarize(across(-contains("weights"), ~first(.)),
              weights = sum(weights),
              subclass = NA,
              dist = NA)
}


agg_avg_units <- function(scweights) {
  if (!is.data.frame(scweights)) {
    scweights <- scweights %>%
      map_dfr(~mutate(., subclass=id[1]))
  }

  scweights %>%
    group_by(subclass, Z) %>%
    summarize(across(starts_with("X"),
                     ~sum(.x * 1/n())),
              across(starts_with("Y"),
                     ~sum(.x * 1/n()))) %>%
    mutate(id = c(NA, subclass[1]), .before="subclass") %>%
    mutate(weights = 1) %>%
    ungroup()
}




