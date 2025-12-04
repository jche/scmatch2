library(latex2exp)


# functions to generate diagnostic plots
create_toy_df_plot <- function(toy_df) {
  toy_df_plot <- toy_df %>%
    mutate(Z = factor(Z, levels = c(FALSE, TRUE), labels = c("Control", "Treatment"))) %>%
    ggplot(aes(x = X1, y = X2, col=as.factor(Z))) +
    geom_point() +
    # scale_color_continuous(low = "blue", high = "orange") +
    theme_classic() +
    labs(x = TeX("$X_1$"),
         y = TeX("$X_2$"),
         color = "Group")
  return(toy_df_plot)
}


# Exploring distances -----


#' Distance of kth neighbor plot
#'
#' Plot distribution of distances between treated units and their
#' kth-nearest neighbors in the control group, for a set of k values.
#'
#' @param scm A csm_matches object
#' @param tops A vector of integers indicating which nearest neighbors
#'   to plot
#' @param caliper Optional caliper value to plot as a vertical line.
#'   Otherwise it will take caliper from scm object.  If NA will plot
#'   no line.
#' @param target_percentile Optional target percentile for caliper,
#'   which will calculate caliper to achieve.
#'
#' @return The plot with some extra attributes.  First is "distances",
#'   the table of distances to the kth nearest neighbor. Second, if
#'   caliper provided, "table" with the table of proportion of
#'   distances below the caliper for each k.  Last is "caliper", the
#'   caliper value used.
#'
#' @export
caliper_distance_plot <- function( scm, tops = 1:3, caliper = NULL,
                                   target_percentile = NULL,
                                   target_k = 1 ) {

  tops = sort(tops)

  # Extract distance matrix from slot
  dist_matrix <- data.frame(t(as.matrix(scm$dm_uncapped)))
  dm_col_sorted <- apply(dist_matrix, 2, sort)

  # Updated plotting function to show the "hard to match" units
  plot_dm <- function(dist_to_plot){
    tibble(d = as.numeric(as.matrix(dist_to_plot))) %>%
      filter(d < 10) %>% # Filter out exact match dummies (1000), keep adaptive ones (~2.6)
      ggplot(aes(d)) +
      geom_histogram(color="black", binwidth=0.1) +
      geom_vline(xintercept = lalonde_params$caliper, col="red") +
      theme_classic() +
      labs(y=NULL, x = TeX("$d(X_t, X_j)$")) +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  }

  # Top 1, 2, 3 distances
  dists <- dm_col_sorted[tops,] %>%
    t() %>%
    as.data.frame() %>%
    set_names( tops ) %>%
    pivot_longer( cols=everything(),
                  names_to = "Rank",
                  values_to = "Distance") %>%
    mutate( Rank = as.numeric( Rank ) )


  plt <- ggplot( dists, aes( Distance )  ) +
    facet_wrap( ~Rank, ncol=1 ) +
    geom_histogram( color="black" ) +
    labs( y = "" ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  if ( !is.null( target_percentile ) ) {
    if( !( target_k %in% tops ) ) {
      stop( "target_k must be one of the values in tops" )
    }

    caliper <- quantile( dists$Distance[ dists$Rank == target_k ],
                         probs = target_percentile )
  }
  attr( plt, "distances" ) <- dists

  if ( is.null( caliper ) ) {
    caliper = params(scm)$caliper
  }

  if ( !is.null( caliper ) && !is.na( caliper) ) {

    tbl <- dists %>%
      group_by( Rank ) %>%
      summarise( pct_below_caliper = mean( Distance <= caliper ) )

    #  ach_percentile <- tbl$pct_below_caliper[ tbl$Rank == target_k ]

    x_loc = max( dists$Distance ) * 0.95

    plt <- plt +
      geom_vline( xintercept = caliper, col="red" ) +
      geom_text(
        data = tbl,
        aes(
          x = x_loc, y = Inf,
          label = scales::percent(pct_below_caliper, accuracy = 1)
        ),
        hjust = 1.1, vjust = 1.1,
        inherit.aes = FALSE
      ) +
      labs( caption = glue::glue( "Caliper = {round( caliper, 3 )}" ) )

    attr( plt, "caliper" ) <- caliper
    attr( plt, "table" ) <- tbl
  }

  plt
}


#' Calculate distances from all matched treatment units to controls
#'
#' Look at distances between each treated unit and their synthetic
#' control, average control, and closest control.
#'
#' @param scm A csm_matches object
#' @param long_table If TRUE, return a long-form table with method and
#'   distance columns.  If FALSE each tx unit is a row.
#'
#' @return A data frame with distances from each treated unit to their
#'   synthetic control, average control, and closest control.
#'
#' @export
get_distance_table <- function( scm,
                                long_table = FALSE ) {
  d = result_table(scm)
  scaling = params(scm)$scaling
  metric = params(scm)$metric
  covariates = params(scm)$covariates
  treatment = params(scm)$treatment
  outcome = params(scm)$outcome

  if ( is.null( scm$treatment_table ) ) {
    scm$treatment_table <- make_treatment_table(scm)
  }

  # check distances bw each tx/sc pair
  sc_dists <- scm %>%
    agg_sc_units()
  stopifnot( all( sc_dists$id[sc_dists[[treatment]]==1] == scm$treatment_table$id) )
  sc_dists <- sc_dists %>%
    gen_dm(scaling = scaling,
           covs = covariates,
           treatment = treatment,
           metric = metric) %>%
    diag()

  avg_dists <- scm %>%
    agg_avg_units()
  stopifnot( all( avg_dists$id[avg_dists[[treatment]]==1] == scm$treatment_table$id) )

  avg_dists <- avg_dists %>%
    gen_dm(scaling = scaling,
           covs = covariates,
           treatment = treatment,
           metric = metric) %>%
    diag()

  nn_dists <- scm %>%
    result_table() %>%
    group_by(Z,subclass) %>%
    filter(dist == min(dist)) %>%
    mutate(weights = 1/n()) %>%
    ungroup() %>%
    agg_avg_units( covariates = covariates,
                   treatment = treatment,
                   outcome = outcome )
  stopifnot( all( nn_dists$id[nn_dists[[treatment]]==1] == scm$treatment_table$id) )
  nn_dists <- nn_dists %>%
    gen_dm(scaling = scaling,
           covs = covariates,
           treatment = treatment,
           metric = metric) %>%
    diag()

  res_list <- list(sc = sc_dists, avg = avg_dists, nn = nn_dists)
  dist_table <- tibble(id = scm$treatment_table$id,
                       SCM = sc_dists,
                       average = avg_dists,
                       closest = nn_dists)

  dist_table <- dist_table %>%
    left_join( scm$treatment_table %>%
                 dplyr::select( id, feasible, matched ),
               by="id" )

  if ( long_table ) {
    dist_table <- dist_table %>%
      pivot_longer( c( `SCM`, `average`, `closest` ),
                  names_to = "method", values_to="distance" ) %>%
    mutate(method = factor( method, levels = c( "SCM", "Average", "1-NN" ) ) )
  }


  dist_table
}


# SCM evaluation plot -----------------------------------------------------

#' Calculate distances from treated units to their controls
#'
#'
#' Look at distribution of pairwise distances between treated unit and
#' their synthetic control, average control, and 1-NN control.
#'
#' @param scm A csm_matches object
#' @param feasible_only If TRUE, only plot distances for treated units
#'   that were feasible and matched.
#' @param boxplot_style If TRUE, use boxplot style for the density
#'   plot.
#'
#' @return A ggplot object showing the density of distances. Also has
#'   two attributes: "table" with summary statistics of the distances,
#'   and "dist_table" with the full table of distances.
#'
#' @export
distance_density_plot <- function(scm, feasible_only = FALSE, boxplot_style = TRUE ) {

  dist_table <- get_distance_table( scm, long_table = TRUE )


  dd <- if ( feasible_only ) {
    dist_table %>%
      filter( feasible == 1 & matched == 1 )
  } else {
    dist_table
  }

  if ( boxplot_style ) {
    plt <- ggplot( dd, aes( method, distance, col=method ) ) +
      geom_boxplot() +
      coord_flip() +
      theme_classic() +
      scale_color_manual(values = wesanderson::wes_palette("Zissou1", 5)[c(5,4,1)]) +
      labs(y = "",
           y = "Distance between treated unit \nand corresponding control",
           color = "Method")

  } else {
    plt <- ggplot(dd, aes(x=distance, color=method)) +
      geom_density(linewidth=1) +
      theme_classic() +
      scale_color_manual(values = wesanderson::wes_palette("Zissou1", 5)[c(5,4,1)]) +
      labs(y = "Density",
           x = "Distance between treated unit \nand corresponding control",
           color = "Method")
  }

  tbl <- dd %>%
    group_by( method ) %>%
    summarise( mean = mean( distance ),
               median = median( distance ) )

  attr( plt, "table" ) <- tbl
  attr( plt, "dist_table" ) <- dist_table

  return(plt)
}



#' Explore pairwise distances
#'
#' Look at the distance between each treated unit and their synthetic control
#' versus the distance between each treated unit and the simple average control.
#'
#' @noRd
scm_vs_avg_distance_plot <- function(scm) {

  d = result_table(scm)
  scaling = params(scm)$scaling
  metric = params(scm)$metric

  # check distances bw each tx/sc pair
  sc_dists <- d %>%
    agg_sc_units() %>%
    gen_dm(scaling = scaling,
           metric = metric) %>%
    diag()
  avg_dists <- d %>%
    agg_avg_units() %>%
    gen_dm(scaling = scaling,
           metric = metric) %>%
    diag()


  # identify places where average does better than sc
  #  (this should never happen)
  if (F) {
    ids <- feasible %>%
      agg_sc_units() %>%
      filter(!is.na(id)) %>%
      pull(id)
    tibble(sc_dists = sc_dists,
           avg_dists = avg_dists) %>%
      mutate(id = ids) %>%
      filter(sc_dists > avg_dists) %>%
      print(n=30)

    foo <- 153
    feasible %>%
      filter(subclass == foo)
    feasible %>%
      agg_sc_units() %>%
      filter(subclass == foo)
    feasible %>%
      agg_avg_units() %>%
      filter(subclass == foo)
  }

  # TODO: simple plot comparing distance between:
  #  - tx unit and simple average control
  #  - tx unit and synthetic control

  tibble(sc_dists = sc_dists,
         avg_dists = avg_dists) %>%
    ggplot(aes(x=sc_dists, y=avg_dists)) +
    geom_point() +
    geom_abline(lty="dotted") +
    theme_classic() +
    labs(y = "Distance between treated and average unit",
         x = "Distance between treated and SC unit")
}


# estimate-estimand tradeoff plot -----------------------------------------

satt_plot <- function(scm, B=NA) {

  feasible_subclasses <- feasible_unit_subclass(scm)

  n_feasible <- length(feasible_subclasses)

  res = result_table(scm)
  ggd_maxcal <- res %>%
    filter(!subclass %in% feasible_subclasses) %>%
    filter(Z==1) %>%
    left_join( caliper_table(scm), by="id") %>%
    arrange(adacal) %>%
    mutate(order = 1:n() + n_feasible)

  p1 <- ggd_maxcal %>%
    ggplot(aes(x=order, y=adacal)) +
    geom_col(fill="gray") +
    geom_hline(yintercept=CALIPER, lty="dotted") +
    theme_classic() +
    labs(y = "Largest adaptive \ncaliper used",
         x = "Total number of treated units used")

  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(caliper_table(scm), by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(cum_avg = cumsum(tx) / (1:n()) ) %>%
    slice((length(attr(res, "feasible_units"))+1):n()) %>%
    mutate(order = 1:n() + n_feasible)

  p2 <- ggd_att %>%
    ggplot(aes(x=order, y=cum_avg)) +
    geom_point() +
    geom_step(linewidth=1) +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Total number of treated units used")

  if (!is.na(B)) {
    # bootstrap dists of FSATT and SATT
    boot_fsatt <- attr(res, "scweights") %>%
      bind_rows() %>%
      filter(subclass %in% feasible_subclasses) %>%
      agg_co_units() %>%
      boot_bayesian(B=B)

    boot_satt <- attr(res, "scweights") %>%
      agg_co_units() %>%
      boot_bayesian(B=B)

    # # using sd of bootstrap samples
    # boot_df <- ggd_att %>%
    #   filter(order == min(order) | order == max(order)) %>%
    #   mutate(sd = c(sd(boot_fsatt), sd(boot_satt)))
    #
    # p2 <- p2 +
    #   geom_errorbar(data=boot_df,
    #                 aes(ymin=cum_avg-2*sd, ymax=cum_avg+2*sd),
    #                 width = 1)

    # using full distribution of bootstrap samples
    boot_df <- ggd_att %>%
      filter(order == min(order) | order == max(order)) %>%
      mutate(q025 = c(quantile(boot_fsatt, 0.025),
                      quantile(boot_satt, 0.025)),
             q975 = c(quantile(boot_fsatt, 0.975),
                      quantile(boot_satt, 0.975)))

    p2 <- p2 +
      geom_errorbar(data=boot_df,
                    aes(ymin=q025, ymax=q975),
                    width = 1)
  }

  p1/p2
}

satt_plot2 <- function(res, B=NA) {
  feasible_subclasses <- attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)

  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(attr(res, "adacalipers"), by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(order = 1:n(),
           cum_avg = cumsum(tx) / order ) %>%
    slice((n_feasible+1):n())

  foo <- ggd_att %>%
    group_by(adacal) %>%
    summarize(cum_avg = last(cum_avg),
              n = n()) %>%
    mutate(shape = ifelse(n==1, 20, 19))
  p <- ggd_att %>%
    ggplot(aes(x=adacal, y=cum_avg)) +
    geom_point(data = foo,
               aes(shape = shape),
               size=2) +
    geom_step() +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Maximum caliper size used") +
    scale_shape_identity()

  if (!is.na(B)) {
    # bootstrap dists of FSATT and SATT
    boot_fsatt <- attr(res, "scweights") %>%
      bind_rows() %>%
      filter(subclass %in% feasible_subclasses) %>%
      agg_co_units() %>%
      boot_bayesian(B=B)

    boot_satt <- attr(res, "scweights") %>%
      agg_co_units() %>%
      boot_bayesian(B=B)

    # # using sd of bootstrap samples
    # boot_df <- ggd_att %>%
    #   filter(order == min(order) | order == max(order)) %>%
    #   mutate(sd = c(sd(boot_fsatt), sd(boot_satt)))
    #
    # p2 <- p2 +
    #   geom_errorbar(data=boot_df,
    #                 aes(ymin=cum_avg-2*sd, ymax=cum_avg+2*sd),
    #                 width = 1)

    # using full distribution of bootstrap samples
    boot_df <- ggd_att %>%
      filter(order == min(order) | order == max(order)) %>%
      mutate(q025 = c(quantile(boot_fsatt, 0.025),
                      quantile(boot_satt, 0.025)),
             q975 = c(quantile(boot_fsatt, 0.975),
                      quantile(boot_satt, 0.975)))

    p <- p +
      geom_errorbar(data=boot_df,
                    aes(ymin=q025, ymax=q975),
                    width = 0.05)
  }

  p
}

satt_plot3 <- function(res, B=NA) {
  feasible_subclasses <- attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)

  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(attr(res, "adacalipers"), by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(order = 1:n(),
           cum_avg = cumsum(tx) / order ) %>%
    slice((n_feasible):n())

  p <- ggd_att %>%
    ggplot(aes(x=order, y=cum_avg)) +
    geom_line(alpha=0.5) +
    geom_point(aes(color=adacal),
               size=3) +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Total number of treated units used",
         color = "Maximum caliper size used") +
    scale_color_continuous(low="blue", high="orange") +
    expand_limits(color=1)

  if (!is.na(B)) {
    # bootstrap dists of FSATT and SATT
    boot_fsatt <- attr(res, "scweights") %>%
      bind_rows() %>%
      filter(subclass %in% feasible_subclasses) %>%
      agg_co_units() %>%
      boot_bayesian(B=B)

    boot_satt <- attr(res, "scweights") %>%
      agg_co_units() %>%
      boot_bayesian(B=B)

    # using sd of bootstrap samples
    boot_df <- ggd_att %>%
      filter(order == min(order) | order == max(order)) %>%
      mutate(sd = c(sd(boot_fsatt), sd(boot_satt)))

    print(paste0("FSATT: (",
                 round(boot_df$cum_avg[1]-1.96*boot_df$sd[1],3),
                 ", ",
                 round(boot_df$cum_avg[1]+1.96*boot_df$sd[1], 3), ")"))
    print(paste0("SATT: (",
                 round(boot_df$cum_avg[2]-1.96*boot_df$sd[2], 3),
                 ", ",
                 round(boot_df$cum_avg[2]+1.96*boot_df$sd[2], 3), ")"))

    p <- p +
      geom_errorbar(data=boot_df,
                    aes(ymin=cum_avg-1.96*sd, ymax=cum_avg+1.96*sd),
                    width = 0.5,
                    linewidth=1) +
      geom_point(aes(color=adacal),
                 size=3)

    # # using full distribution of bootstrap samples
    # boot_df <- ggd_att %>%
    #   filter(order == min(order) | order == max(order)) %>%
    #   mutate(q025 = c(quantile(boot_fsatt, 0.025),
    #                   quantile(boot_satt, 0.025)),
    #          q975 = c(quantile(boot_fsatt, 0.975),
    #                   quantile(boot_satt, 0.975)))
    #
    # p <- p +
    #   geom_errorbar(data=boot_df,
    #                 aes(ymin=q025, ymax=q975),
    #                 width = 0.5) +
    #   geom_point(aes(color=adacal),
    #              size=3)
  }

  p
}

# fixed, with wild bootstrap
satt_plot4 <- function(res, B=NA) {
  feasible_subclasses <-
    attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)

  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(attr(res, "adacalipers"), by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(order = 1:n(),
           cum_avg = cumsum(tx) / order ) %>%
    slice((n_feasible):n())

  p <- ggd_att %>%
    ggplot(aes(x=order, y=cum_avg)) +
    geom_line(alpha=0.5) +
    geom_point(
      # aes(color=adacal),
      size=3) +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Total number of treated units used",
         color = "Maximum caliper size used") +
    # scale_color_continuous(low="blue", high="orange") +
    expand_limits(color=1)

  if (!is.na(B)) {
    # bootstrap dists of FSATT and SATT
    boot_fsatt <- res %>%
      filter(subclass %in% feasible_subclasses) %>%
      wild_bootstrap(B=B)

    boot_satt <- attr(res, "scweights") %>%
      agg_sc_units() %>%
      wild_bootstrap(B=B)

    # using sd of bootstrap samples
    boot_df <- ggd_att %>%
      filter(order == min(order) | order == max(order)) %>%
      mutate(sd = c(sd(boot_fsatt), sd(boot_satt)))

    print(paste0("FSATT: (",
                 round(boot_df$cum_avg[1]-1.96*boot_df$sd[1],3),
                 ", ",
                 round(boot_df$cum_avg[1]+1.96*boot_df$sd[1], 3), ")"))
    print(paste0("SATT: (",
                 round(boot_df$cum_avg[2]-1.96*boot_df$sd[2], 3),
                 ", ",
                 round(boot_df$cum_avg[2]+1.96*boot_df$sd[2], 3), ")"))

    p <- p +
      geom_errorbar(data=boot_df,
                    aes(ymin=cum_avg-1.96*sd, ymax=cum_avg+1.96*sd),
                    width = 0.5,
                    linewidth=1) +
      geom_point(aes(color=adacal),
                 size=3)
  }

  p
}


# love plot ---------------------------------------------------------------

get_diff_scm_co_and_tx <- function(res, covs){
  stopifnot( is.csm_matches(res) )

  ada = res$treatment_table %>%
    select(id, adacal) %>%
    mutate( id = as.character(id) )
  df_diff_scm_co_and_tx <- result_table(res) %>%
    left_join(ada,
              by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              across(all_of(covs),
                     ~.[2] - .[1]))
  return(df_diff_scm_co_and_tx)
}



create_love_plot_df <- function(res, covs){
  feasible_subclasses <- feasible_units(res)
  n_feasible <- nrow(feasible_subclasses)

  df_step_1<-
    get_diff_scm_co_and_tx(res,covs)

  df_step_2<-df_step_1 %>%
    arrange(adacal) %>%
    mutate(order = 1:n(),
           across(all_of(covs),
                  ~cumsum(.) / order))

  df_step_3 <- df_step_2 %>%
    slice((n_feasible):n()) %>%
    pivot_longer(all_of(covs))
  love_plot_df <- df_step_3
  return(love_plot_df)
}


#' Love plot of covariate balance
#'
#' Make a ggplot love plot of covariate balance for each covariate
#' passed.
#'
#' @export
love_plot <- function( csm, covs = NULL, covs_names = NULL, B=NA ) {


  cc = params(csm)$covariates
  if ( !is.null( covs ) ) {
    stopifnot( all( covs %in% cc ) )
  } else {
    covs = cc
  }

  if ( !is.null( covs_names ) && !is.null( covs ) ) {
    stopifnot( length( covs_names ) == length( covs ) )
  }

  love_steps <- create_love_plot_df(csm, covs)

  if ( !is.null( covs_names ) ) {
    love_steps$name <- covs_names[ match( love_steps$name, covs ) ]
  } else {
    covs_names = covs
  }

  p <- love_steps %>%
    ggplot(aes(x=order, color=name)) +
    geom_point(data=. %>%
                 slice( c( seq(1,length(covs)),
                           seq( n()-length(covs)+1, n()))),
               aes(y=value), size=2) +
    geom_step(aes(y=value, group=name),
              linewidth=1.1) +
    expand_limits(y = 0) +
    geom_hline(yintercept=0,
               lty="dotted") +
    facet_wrap(~name,
               scales="free_y") +
    labs(y = "\n Covariate balance (tx-co)",
         x = "Total number of treated units used",
         color = "Covariate") +
    theme_classic() +
    make_tx_axis( love_steps$order )

    if (!is.na(B)) {
      # bootstrap dists of FSATT and SATT
      boot_fsatt <- attr(csm, "scweights") %>%
        bind_rows() %>%
        filter(subclass %in% feasible_subclasses) %>%
        agg_co_units() %>%
        boot_bayesian_covs(covs=covs, B=B) %>%
        pivot_longer(everything()) %>%
        group_by(name) %>%
        summarize(q025 = quantile(value, 0.025),
                  q975 = quantile(value, 0.975)) %>%
        mutate(order = n_feasible+1)

      boot_satt <- attr(csm, "scweights") %>%
        agg_co_units() %>%
        boot_bayesian_covs(covs=covs, B=B) %>%
        pivot_longer(everything()) %>%
        group_by(name) %>%
        summarize(q025 = quantile(value, 0.025),
                  q975 = quantile(value, 0.975)) %>%
        mutate(order = max(love_steps$order))

      p <- p +
        geom_errorbar(data=bind_rows(boot_fsatt, boot_satt),
                      aes(ymin=q025, ymax=q975),
                      width = 1)
    }

  return(p)

}


love_plot2 <- function(res, covs, B=NA) {
  feasible_subclasses <- attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)

  adacal_key <- attr(res, "adacalipers") %>%
    rename(subclass = id) %>%
    arrange(adacal) %>%
    mutate(order = 1:n())

  love_steps <- res %>%
    left_join(adacal_key, by="subclass") %>%
    group_by(Z) %>%
    arrange(order) %>%
    mutate(across(all_of(covs), ~cumsum(.) / order)) %>%
    slice((n_feasible):n()) %>%
    pivot_longer(all_of(covs))

  love_steps %>%
    mutate(Z = as.factor(Z)) %>%
    ggplot(aes(x=order, y=value, color=Z)) +
    # geom_line(aes(group=order), alpha=0.3, color="black") +
    geom_point(size=1) +
    geom_point(data=. %>%
                 filter(order == max(order) |
                          order == min(order)),
               size=2) +
    geom_line(aes(group=interaction(name, Z)),
              alpha=0.3) +
    facet_wrap(~name, scales="free_y") +
    labs(y = "Marginal mean",
         x = "Total number of treated units used",
         color = "Z") +
    theme_classic()
}

