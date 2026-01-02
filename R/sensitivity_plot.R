



#' Make sensitivity table of different ATT estimates
#'
#' Given a csm_matches object, this function computes several different
#' ATT estimates:
#' - FATT: ATT using only feasible treated units (i.e., those with
#'  matches within the caliper)
#'  - ATT: ATT using all treated units, if there were adaptive calipers.
#'  - ATT_1nn: ATT using nearest-neighbor controls (if tied, use all tied controls).
#'  - ATT_raw: ATT, averaging all controls within each matched set.
#'  - ATT_raw_feas: ATT, averaging all controls within each matched set, but only for feasible treated units.
#'
#' These estimates can be used to assess sensitivity of results to
#' different choices in matching and estimation.  Return results are from the `estimate_ATT()` method
#'
#'
#' @seealso [estimate_ATT()]
#'
#' @param csm A csm_matches object
#' @param outcome Name of outcome variable in data
#' @return A tibble with the four ATT estimates and standard errors
#'
#' @export
sensitivity_table <- function( csm,
                               outcome = NULL,
                               feasible_only = FALSE,
                               include_variances = FALSE,
                               include_distances = TRUE ) {

  args = attr( csm, "settings" )


  if ( !feasible_only && args$rad_method != "adaptive" ) {
    stop( "Before calling sensitivity_table, refit csm object to adaptive radii so sensitivity_table can obtain full result table" )
  }
  if ( feasible_only && !( args$rad_method %in% c( "fixed", "adaptive" ) ) ) {
    stop( "Before calling sensitivity_table, refit csm object to adaptive or fixed radii so sensitivity_table can obtain full result table" )
  }

  tx_var = args$treatment
  covs = attr( csm, "covariates" )

  if ( args$est_method != "scm" ) {
    csm <- est_weights( matched_gps = csm,
                        est_method = "scm")
  }

  FATT = estimate_ATT( csm,
                       treatment = tx_var,
                       outcome = outcome,
                       feasible_only = TRUE )


  # 1nn estimate
  rs = result_table( csm, include_caliper = TRUE )

  rs_knn <- rs %>%
    group_by( subclass, .data[[tx_var]] ) %>%
    filter( dist == min(dist) ) %>%
    mutate( weights = 1 / n() ) %>%
    ungroup()
  ATT_1nn_feas = estimate_ATT( rs_knn %>% filter( feasible == 1 ),
                               treatment = tx_var,
                               outcome = outcome )

  rs_raw <- rs %>%
    group_by( subclass, .data[[tx_var]] ) %>%
    mutate( weights = 1 / n() ) %>%
    ungroup()
  ATT_raw_feas = estimate_ATT( rs_raw %>% filter( feasible == 1 ),
                               treatment = tx_var,
                               outcome = outcome )

  if ( !feasible_only ) {

    ATT = estimate_ATT( csm,
                        treatment = tx_var,
                        outcome = outcome,
                        feasible_only = FALSE )
    ATT_1nn = estimate_ATT( rs_knn,
                            treatment = tx_var,
                            outcome = outcome )
    ATT_raw = estimate_ATT( rs_raw,
                            treatment = tx_var,
                            outcome = outcome )

    rs <- bind_rows( FATT = FATT,
                     ATT = ATT,
                     FATT_1nn = ATT_1nn_feas,
                     ATT_1nn = ATT_1nn,
                     FATT_raw = ATT_raw_feas,
                     ATT_raw = ATT_raw,
                     .id = "Estimate" )
  } else {
    rs <- bind_rows( FATT = FATT,
                     FATT_1nn = ATT_1nn_feas,
                     FATT_raw = ATT_raw_feas,
                     .id = "Estimate" )
  }

  if ( !is.null(outcome) && !include_variances ) {
    rs <- rs %>%
      dplyr::select( -V, -V_E, -V_P )
  }

  if ( include_distances ) {
    ddp = distance_density_plot( csm, feasible_only = FALSE )
    dplt <- attr( ddp, "table" )
    dplt_feas = attr( ddp, "dist_table" ) %>%
      dplyr::filter( feasible == 1 ) %>%
      group_by( method ) %>%
      summarise( mean = mean( distance ),
                 median = median( distance ) )

    all_dists = bind_rows( dplt, dplt_feas )
    all_dists$method = c( "ATT", "ATT_raw", "ATT_1nn",
                          "FATT", "FATT_raw", "FATT_1nn" )
    all_dists <- all_dists %>%
      rename( mean_dist = mean,
              median_dist = median )

    rs <- rs %>%
      left_join( all_dists, by = c( "Estimate" = "method" ) ) %>%
      select( Estimate, everything() )

  }

  rs <- arrange( rs, Estimate )

  if ( feasible_only ) {
    target = "FATT_raw"
  } else {
    target = "ATT_raw"
  }

  rs <- rs %>%
    mutate(
      SE_star = sqrt( 1/N_T + 1/ESS_C )
    ) %>%
    mutate( SE_ratio = SE_star / SE_star[ Estimate == target ],
            bias_ratio = mean_dist / mean_dist[ Estimate == target ] )

  #    ratio = mean_dist / SE_star ) %>%
  #  mutate( ratio = ratio / ratio[ Estimate == target ] )

  rs
}




#' Make sensitivity table or plot of impact of changing caliper
#'
#' To generate the sensitivity table, method will run caliper from a
#' min to max value, generating matches at each step. Then make a plot
#' of the ATT and confidence interval for each point on the plot. Also
#' can return table of results if 'return_table' is set to TRUE.
#'
#'
#' @rdname caliper_sensitivity_table
#'
#' @details `caliper_sensitivity_table()` generates the sensitivity
#'   table only.
#'
#' @param csm A CSM object OR a sensitivity ggplot created by, e.g.,
#'   [caliper_sensitivity_plot()] or
#'   [caliper_sensitivity_plot_stats()] (or anything with the
#'   sensitivity table as an attribute) OR a table with the relevant
#'   columns (focus, lines, caliper, ESS_C, etc.). This allows calls
#'   to the different plot methods without recalculating everything.
#' @param data The data the csm_matches object was fit to. (One could
#'   put in alternate data here, if the covariates aligned.)
#' @param min_cal Minimum caliper to try
#' @param max_cal Maximum caliper to try
#' @param R Number of calipers to try between min and max
#'
#' @return A tibble of results.
#'
#' @export
caliper_sensitivity_table <- function( csm,
                                       data,
                                       outcome = NULL,
                                       include_distances = TRUE,
                                       min_cal = 0.05,
                                       max_cal = 5,
                                       R = 30 ) {

  # Short circuit if the object happens to be prior plot or something
  # with a table attribute.
  if ( is.data.frame( csm ) ) {
    #stopifnot( all( c( "caliper", "Estimate", "ATT", "SE", "ESS_C" ) %in% colnames(csm) ) )
    return( csm )
  }
  st = attr( csm, "table" )
  if ( !is.null(st) ) {
    return( st )
  }

  cals = seq( min_cal, max_cal, length.out = R )

  args = attr( csm, "settings" )
  covs = attr( csm, "covariates" )

  tx_var = args$treatment

  if ( params(csm)$rad_method != "adaptive"  ) {
    #warn( "Refitting csm object to adaptive radii to calculate full range of estimates" )
    csm <- update_matches( csm, data,
                           rad_method = "adaptive" )
  }


  sen <- map( cals, function( cali ) {
    update_matches( csm, data=data,
                    caliper = cali )
  } )


  sens_table <- sen %>%
    set_names( cals ) %>%
    purrr::map_dfr( sensitivity_table,
                    outcome = outcome,
                    include_distances = include_distances,
                    .id = "caliper" ) %>%
    dplyr::mutate( caliper = as.numeric( caliper ) ) %>%
    dplyr::relocate( caliper )


  sens_table

}




pick_lines <- function( sens_table ) {

  sens_table <- sens_table %>%
    mutate( feasible = ifelse( Estimate %in% c("FATT","FATT_raw", "FATT_1nn"),
                               "Feasible", "All" ),
            weight = dplyr::case_when(
              grepl("_1nn$", Estimate) ~ "1nn",
              grepl("_raw$", Estimate) ~ "Avg",
              TRUE                 ~ "Synth"
            ) )

  sens_table
}


sens_table_legend <- function() {
  list( scale_color_manual( values = c( "Feasible" = "blue",
                                        "All" = "orange" ) ),
        scale_linetype_manual( values = c( "Synth" = "solid",
                                           "Avg" = "dashed",
                                           "1nn" = "dotted" ) ),
        theme_minimal(),
        labs( y = "Estimate", x = "Caliper",
              color = "Group",
              linetype = "Weight" )
  )
}


#' @rdname caliper_sensitivity_table
#'
#' @details `caliper_sensitivity_plot()` generates a ggplot with
#' y-axis as ATT estimate and confidence interval, and x-axis a series
#' of calipers from the minimum to maximum specified.  Different lines
#' are shown for different ATT estimates.
#'
#' @param focus What type of estimate to plot the confidence envelope
#'   for.
#' @param lines Which estimates to show as lines on the plot.
#'
#' @export
caliper_sensitivity_plot <- function( csm,
                                      data,
                                      outcome = NULL,
                                      include_distances = TRUE,
                                      min_cal = 0.05,
                                      max_cal = 5,
                                      R = 30,
                                      focus = NULL,
                                      lines = c( "FATT", "ATT", "FATT_1nn", "ATT_1nn",
                                                 "FATT_raw", "ATT_raw" ) ) {

  sens_table <- caliper_sensitivity_table( csm,
                                           data,
                                           outcome = outcome,
                                           include_distances = include_distances,
                                           min_cal = min_cal,
                                           max_cal = max_cal,
                                           R = R )

  if ( !( "ATT" %in% colnames(sens_table) ) ) {
    if ( is.csm_matches( csm ) ) {
      stop( "ATT estimate not found in csm_matches object---rerun sensitivity table with specified outcome" )
    } else {
      stop( "ATT estimate not found in table---rerun sensitivity table with specified outcome" )
    }
  }

  stopifnot( all( c( "caliper", "Estimate", "ATT", "SE", "ESS_C" ) %in% colnames(sens_table) ) )


  sens_table = pick_lines( sens_table )
  #focus = match.arg( focus )

  lines = filter( sens_table, Estimate %in% lines )

  plt <- ggplot( lines, aes( x=caliper,
                             y=ATT ) ) +
    geom_hline( yintercept = 0, color="black" )

  if ( !is.null( focus ) ) {
    atts = filter( sens_table, Estimate %in% focus )

    plt <- plt +
      geom_ribbon( data=atts, aes( ymin=ATT-2*SE, ymax=ATT+2*SE ), alpha=0.2 ) +
      geom_line( data=atts,
                 aes( col = feasible,
                      linetype = weight,
                      group = Estimate ),
                 lwd=2 )
  }

  plt <- plt + geom_line( data=lines,
                          aes( col = feasible,
                               linetype = weight,
                               group = Estimate ) ) +
    sens_table_legend()

  attr( plt, "table" ) <- sens_table

  return( plt )

}


#' @rdname caliper_sensitivity_table
#'
#' @details `caliper_sensitivity_plot_stats()` makes an augmented plot
#' from sensitivity plot showing ATT, SE and ESS
#'
#' @inheritParams caliper_sensitivity_plot
#'
#' @export
caliper_sensitivity_plot_stats <- function( csm,
                                            data,
                                            outcome = NULL,
                                            min_cal = 0.05,
                                            max_cal = 5,
                                            R = 30,
                                            vars = c( "ESS_C", "mean_dist", "SE_star" ),
                                            focus = NULL,
                                            lines = c( "FATT", "ATT", "FATT_1nn", "ATT_1nn",
                                                       "FATT_raw", "ATT_raw" ) ) {

  sens_table <- caliper_sensitivity_table( csm,
                                           data,
                                           outcome = outcome,
                                           min_cal = min_cal,
                                           max_cal = max_cal,
                                           R = R )


  sens_table = pick_lines( sens_table )

  tblL <- sens_table %>%
    pivot_longer( cols = all_of( vars ),
                  names_to = "Metric",
                  values_to = "Value" )

  tblL$Metric = factor( tblL$Metric, levels = vars,
                        ordered = TRUE )

  lines = filter( tblL, Estimate %in% lines )

  cc = sens_table$caliper[[1]]

  # Comparing ATT, SE and effective sample sizes
  plt <- ggplot( lines, aes( caliper, Value,
                             col = feasible, linetype=weight, group = Estimate ) ) +
    facet_wrap( ~ Metric, scales = "free_y", nrow=1 ) +
    geom_line() +
    #geom_blank(data = data.frame(Metric = "ESS_C", feasible = "All", Estimate = "FATT", weight="Synth", caliper = cc, Value = 0)) +
    #geom_blank(data = data.frame(Metric = "SE", feasible = "All", Estimate = "FATT", weight="Synth", caliper = cc, Value = 0)) +
    geom_hline( yintercept = 0, color="black" ) +
    sens_table_legend()

  if ( !is.null( focus ) ) {
    atts = filter( tblL, Estimate == focus )
    plt <- plt +
      geom_line( data=atts, lwd=2 )
  }

  attr( plt, "table" ) <- sens_table

  return( plt )
}




make_tx_axis <- function( ntxes ) {
  nxtes = sort( unique( ntxes ) )

  dp = length(ntxes)
  rest_pts = ntxes[-seq(1,ceiling(dp/10))]

  mn <- min(ntxes)
  brk <- c( mn, pretty(rest_pts))
  brk = sort( brk )
  lbl = brk
  if ( mn != brk[[1]] ) {
    lbl[1] = ""
  }
  scale_x_continuous(breaks = brk, labels = lbl)

}

#' Make feasible plot showing cumulative ATT as feasible units are
#' added
#'
#' Given a csm_matches object, this function makes a plot showing how
#' the cumulative ATT estimate changes as more feasible treated units
#' are added (i.e., as the caliper size increases to include more
#' treated units).
#'
#' @param csm A csm_matches object
#' @param return_table If TRUE, return the full table of results
#'   instead of a plot
#' @param caliper_plot If TRUE, return a plot of maximum caliper size
#'   vs number of treated units.  If "both" return a list of both
#'   plots, along with the table itself.
#' @return A ggplot object showing the feasible plot, or a table of
#'   results if return_table is TRUE.
#' @export
feasible_plot <- function( csm,
                           outcome = "Y",
                           return_table = FALSE,
                           caliper_plot = FALSE ) {

  tx_var = attr( csm, "settings" )$treatment
  full_table <-
    left_join( result_table( csm ),
               caliper_table( csm ) %>%
                 dplyr::select( -id ),
               by = "subclass" ) %>%
    arrange( max_dist, subclass, -.data[[tx_var]] ) %>%
    dplyr::select( id, subclass, {{outcome}}, {{tx_var}}, weights, feasible, adacal, max_dist )

  full_table <- full_table %>%
    mutate( order = factor( subclass,
                            levels = unique( subclass ),
                            ordered=TRUE ) )

  feasible_subclasses <- feasible_unit_subclass( csm )

  # Make all the feasible units one grouping
  full_table$order[ full_table$feasible == 1 ] = full_table$order[[1]]
  full_table$order = droplevels( full_table$order )

  # Make a table of the att
  ggd_att <- map_dfr( levels( full_table$order ),
                      function( sc ) {
                        df_subclass <- full_table %>%
                          filter( order <= sc )

                        aa <- estimate_ATT( df_subclass,
                                            treatment=tx_var, outcome=outcome )
                        aa$adacal = max( df_subclass$adacal )
                        aa
                      } )

  if ( return_table ) {
    return( ggd_att )
  }

  plot_max_caliper_size <- NA
  if ( caliper_plot == TRUE || (caliper_plot=="both") ) {
    plot_max_caliper_size <-
      ggd_att %>%
      ggplot(aes(x=N_T, y=adacal)) +
      geom_line(alpha=0.5) +
      geom_point(size=3) +
      theme_classic() +
      labs(y = "Maximum caliper size used",
           x = "Total number of treated units used") +
      expand_limits(y=0) +
      make_tx_axis( ggd_att$N_T )

    if ( caliper_plot != "both" ) {
      return( plot_max_caliper_size )
    }
  }



  # PLot the results


  plot_SATT <-  ggd_att %>%
    ggplot(aes(x=N_T, y=ATT)) +
    geom_line(alpha=0.5) +
    geom_point(size=2) +
    geom_point( data=ggd_att[1,], size=5 ) +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Total number of treated units used") +
    #expand_limits(color=1) +
    # theme(
    #   legend.direction="horizontal",
    #    legend.position.inside = c(0.5, 0.85),
    #   legend.background = element_blank(),
    #    legend.box.background = element_rect(colour = "black")
    #  ) +
    labs(
      y = "Cumulative ATT Estimate",
      x = "Total number of treated units used"
    ) +
    geom_ribbon(
      aes(
        ymin = ATT - 1.96 * SE,
        ymax = ATT + 1.96 * SE
      ),
      alpha=0.2
      #width = 0.1,
      #linewidth = 0.5
    ) +
    geom_hline(yintercept=0, lty="dotted") +
    geom_hline(yintercept=ggd_att$ATT[1], lty="dashed", color="blue") +
    make_tx_axis( ggd_att$N_T )

  # ylim(c(-0.1, 0.2))

  if ( caliper_plot == "both" ) {
    return( list( plot_SATT = plot_SATT,
                  plot_max_caliper_size = plot_max_caliper_size,
                  table = ggd_att ) )
  }
  plot_SATT

}


