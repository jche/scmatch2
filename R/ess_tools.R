
# Effective Sample Size exploration tools




# ESS plot ----------------------------------------------------------------

#' Effective sample size (ESS) plot
#'
#' Plot the effective sample size (ESS) of the treated units, and
#' various forms of possible control groups depending on their
#' weighting.
#'
#' @param d Data frame of matched units including columns for Z
#'   (treatment indicator), subclass (matched set identifier), dist
#'   (distance from treated unit), and weights (SCM weights).  E.g., the results from the
#'
#' @export
ess_plot <- function(d) {
  df_sc <- d %>%
    filter(Z==FALSE)

  df_avg <- d %>%
    filter(Z==F) %>%
    group_by(subclass) %>%
    mutate(weights = 1/n()) %>%
    ungroup()

  df_1nn <- d %>%
    filter(Z==F) %>%
    group_by(subclass) %>%
    filter(dist == min(dist)) %>%
    slice(1) %>%
    mutate(weights = 1) %>%
    ungroup()

  df_tx <- d %>%
    filter(Z==T)

  vars_df <- list(df_tx, df_sc, df_avg, df_1nn) %>%
    map_dbl(~ .x %>%
              agg_co_units() %>%
              summarize(ess = sum(weights)^2 / sum(weights^2)) %>%
              pull(ess)) %>%
    tibble(method = c("tx", "sc", "avg", "1nn"),
           ess = .)

  print(vars_df)

  vars_df %>%
    mutate(method = factor(method, levels=c("1nn", "avg", "sc", "tx"))) %>%
    ggplot(aes(x=method, y=ess)) +
    geom_col(aes(fill = method)) +
    scale_fill_manual(values = c("black",
                                 wesanderson::wes_palette("Zissou1", 5)[c(5,4,1)])) +
    scale_x_discrete(labels = c("Control units \n(1-NN weights)",
                                "Control units \n(Average weights)",
                                "Control units \n(SCM weights)",
                                "Treated units")) +
    theme_classic() +
    guides(fill="none") +
    labs(y = "ESS of set of units",
         x = "") +
    coord_flip()
}






#' Compare Effective Sample Size (ESS) across methods
#'
#' Given a csm_matches object, this function computes and compares
#' the Effective Sample Size (ESS) for three different estimation
#' methods: Synthetic Control Method (SCM), Average, and 1-Nearest
#' Neighbor (1-NN).
#'
#' @param csm A csm_matches object
#' @param data The data the csm_matches object was fit to.
#'
#' @return A tibble comparing the ESS for each method
#' @export
compare_ess <- function( csm, data ) {

  # A. SCM (Already run)
  ess_scm <- result_table(csm, feasible_only=TRUE) %>%
    filter(Z==0) %>%
    summarise(ess = sum(weights)^2/sum(weights^2)) %>%
    pull(ess)

  # B. Average (Simple average of matches in caliper)
  csm_avg <-  update_matches( csm,
                              data = data,
                              est_method = "average" )

  ess_avg <- result_table(csm_avg, feasible_only=TRUE) %>%
    filter(Z==0) %>%
    summarise(ess = sum(weights)^2/sum(weights^2)) %>%
    pull(ess)

  # C. 1-NN (Nearest Neighbor)
  # We use rad_method="1nn" and average weighting (since it's 1 unit, avg=1)
  csm_1nn <- update_matches( csm,
                             data = data,
                             rad_method = "1nn",
                             est_method = "average" )

  ess_1nn <- result_table(csm_1nn, feasible_only=TRUE) %>%
    filter(Z==0) %>%
    summarise(ess = sum(weights)^2/sum(weights^2)) %>%
    pull(ess)

  # Print to console to verify against paper (Expect: ~72.6, ~129.9, ~67.0)
  tibble(Method = c("SCM", "Average", "1-NN"),
         ESS = c(ess_scm, ess_avg, ess_1nn))

}



