
# aggregation functions ---------------------------------------------------


#' aggregation methods for matched data
#'
#' agg_sc_units weights to make pairs of tx and synthetic control
#' units. Generally aggregate by cluster defined by treated unit,
#' calculating the weighted average of the control units in the
#' cluster for all variables that start with "X" or "Y" in the
#' dataframes.
#'
#' agg_co_units Aggregates across the controls, grouping repeated
#' controls used by multiple treated units
#'
#' agg_avg_units does simple average of controls within groups (e.g.,
#' like CEM would do).
#'
#' @name AggregationMethods
#' @aliases agg_sc_units, agg_co_units, agg_avg_units
#'
#' @param scweights scm_matches object or list of sc weights (tibbles
#'   with tx unit in first row).
#'
#' @return data.frame with tx and corresponding control units.
#'
#' @export
agg_sc_units <- function(scweights) {
  if ( is.csm_matches(scweights) ) {
    scweights <- scweights$matches
  }

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


#' @rdname AggregationMethods
#' @export
agg_co_units <- function(scweights) {

  if ( is.csm_matches(scweights) ) {
    scweights <- scweights$matches
  }

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



#' @rdname AggregationMethods
#' @export
agg_avg_units <- function(scweights) {
  if ( is.csm_matches(scweights) ) {
    scweights <- scweights$matches
  }

  if (!is.data.frame(scweights)) {
    scweights <- scweights %>%
      map_dfr(~mutate(., subclass=id[1]))
  }

  rs <- scweights %>%
    group_by(subclass, Z) %>%
    summarize(across(starts_with("X"),
                     ~sum(.x * 1/n())),
              across(starts_with("Y"),
                     ~sum(.x * 1/n()))) %>%
    mutate(id = c(NA, subclass[1]), .before="subclass") %>%
    mutate(weights = 1) %>%
    ungroup()

  # Make IDs for the synthetic controls
  rs$id = as.character(rs$id)
  rs$id[ is.na(rs$id) ] <- paste0( rs$subclass[ is.na(rs$id) ], "_avg" )

  rs
}




