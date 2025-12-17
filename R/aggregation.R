
# aggregation functions ---------------------------------------------------


prep_data <- function(scweights,
                      covariates = get_x_vars(scweights),
                      treatment = "Z",
                      outcome =  NULL) {

  # Pull params from the scweights object if it is a csm_matches
  # object
  if ( is.csm_matches(scweights) ) {
    treatment = params(scweights)$treatment
    covariates = attr( scweights, "covariates" )
    scweights <- scweights$matches
  }

  # Convert list of dataframes to single dataframe
  if (!is.data.frame(scweights)) {
    scweights <- scweights %>%
      map_dfr(~mutate(., subclass=id[1]))
  }

  list( scweights = scweights,
        covariates = covariates,
        treatment = treatment,
        outcome = outcome )
}



#' Aggregation methods for matched data
#'
#' @description These functions aggregate matched or
#' synthetic-control–style output into unit-level summaries.
#'
#' @details
#' **`agg_sc_units()`**
#' Aggregate weights within treated-unit subclass to make pairs of tx
#' and corresponding synthetic control units. Generally aggregate by
#' cluster defined by treated unit, calculating the weighted average
#' of the control units in the cluster for all covariates any any
#' provided outcomes
#'
#' **`agg_co_units()`**
#' Aggregates across control units, collecting repeated controls used
#' across treated units and summing their weights.

#'
#' **`agg_avg_units()`**
#' Computes simple averages of controls within each subclass
#' (CEM-style averaging).
#'
#' @param scweights Either a `csm_matches` object or a list/data frame
#'   of matching/synthetic-control weights where each subclass
#'   contains one treated row and its controls.
#' @param covariates Character vector of covariate names. Defaults to
#'   `get_x_vars(scweights)`.
#' @param treatment Name of the treatment indicator column.
#' @param outcome Name of the outcome column.
#'
#' @return A `data.frame` with aggregated treated and control
#'   summaries.
#'
#' @name AggregationMethods
#' @aliases agg_sc_units agg_co_units agg_avg_units
#'
NULL



#' @rdname AggregationMethods
#' @export
agg_sc_units <- function(scweights,
                         covariates = get_x_vars(scweights),
                         treatment = "Z",
                         outcome = NULL ) {
  prep <- prep_data( scweights,
                     covariates,
                     treatment,
                     outcome )
  treatment <- prep$treatment
  outcome <- prep$outcome
  covariates <- prep$covariates
  scweights <- prep$scweights

  cc = c( covariates, outcome )
  rs <- scweights %>%
    group_by(subclass, .data[[treatment]] ) %>%
    summarize(across( all_of( cc ),
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
agg_co_units <- function(scweights,
                         covariates = get_x_vars(scweights),
                         treatment = "Z",
                         outcome = NULL ) {

  prep <- prep_data( scweights,
                     covariates,
                     treatment,
                     outcome )
  treatment <- prep$treatment
  outcome <- prep$outcome
  covariates <- prep$covariates
  scweights <- prep$scweights

  scweights %>%
    group_by(id) %>%
    summarize(across(-contains("weights"), ~first(.)),
              weights = sum(weights),
              subclass = NA,
              dist = NA)
}



#' @rdname AggregationMethods
#' @export
agg_avg_units <- function(scweights,
                          covariates = get_x_vars(scweights),
                          treatment = "Z",
                          outcome = NULL ) {
  prep <- prep_data( scweights,
                     covariates,
                     treatment,
                     outcome )
  treatment <- prep$treatment
  outcome <- prep$outcome
  covariates <- prep$covariates
  scweights <- prep$scweights

  cc = c( covariates, outcome )

  rs <- scweights %>%
    group_by(subclass, Z) %>%
    summarize( across( all_of( cc ),
                     ~sum(.x * 1/n())),
               .groups = "drop_last" ) %>%
    mutate(id = c(NA, subclass[1]), .before="subclass") %>%
    mutate(weights = 1) %>%
    ungroup()

  # Make IDs for the synthetic controls
  rs$id = as.character(rs$id)
  rs$id[ is.na(rs$id) ] <- paste0( rs$subclass[ is.na(rs$id) ], "_avg" )

  rs
}




