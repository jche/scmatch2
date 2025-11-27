
#' Make sensitivity plot of changing caliper
#'
#' Run caliper from a min to max value, generating matches at each
#' step. Then make a plot of the ATT and confidence interval for each
#' point on the plot. Also can return table of results if
#' 'return_table' is set to TRUE.
#'
#' @param res A csm_matches object, used to get the settings for
#'   matching
#' @param data The data the csm_matches object was fit to. (One could
#'   put in alternate data here, if the covariates aligned.)
#' @param min_cal Minimum caliper to try
#' @param max_cal Maximum caliper to try
#' @param R Number of calipers to try between min and max
#' @param return_table If TRUE, return a table of results instead of a
#'   plot
#' @return A ggplot object showing the sensitivity plot.  It has an
#'   attribute of "table" that holds the data used to make the plot.
#' @export
sensitivity_plot <- function( csm,
                              data,
                              min_cal = 0.05,
                              max_cal = 5,
                              R = 30,
                              return_table = FALSE ) {

  cals = seq( min_cal, max_cal, length.out = R )

  args = attr( csm, "settings" )
  covs = attr( csm, "covariates" )

  sen <- map( cals, function( cali ) {
    get_cal_matches( df=data,
                     covs = covs,
                     treatment = args$treatment,
                     metric = args$metric,
                     caliper = cali,
                     rad_method = args$rad_method,
                     est_method = args$est_method,
                     scaling = args$scaling,
                     id_name = args$id_name,
                     k = args$k,
                     warn = FALSE )
  } )

  #sen[[1]]

  atts <- sen %>%
    map_dfr( estimate_ATT ) %>%
    mutate( caliper = cals ) %>%
    relocate( caliper )
  atts

  if ( return_table ) {
    return( atts )
  } else {
    plt <- ggplot( atts, aes( caliper, ATT, ymin=ATT-2*SE, ymax=ATT+2*SE ) ) +
      geom_hline( yintercept = 0, linetype="dashed", color="red" ) +
      geom_ribbon( alpha=0.2 ) +
      geom_line( aes( y=ATT ) ) +
      theme_minimal()

    attr( plt, "table" ) <- atts
    return( plt )
  }
}



make_tx_axis <- function( atts ) {

  dp = nrow(atts)
  rest_pts = atts$N_T[-seq(1,ceiling(dp/10))]

  scale_x_continuous(breaks = c( min(atts$N_T), pretty(rest_pts)))

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
                           return_table = FALSE,
                           caliper_plot = FALSE ) {

  full_table <-
    left_join( result_table( csm ),
               caliper_table( csm ) %>%
                 dplyr::select( -id ),
               by = "subclass" ) %>%
    arrange( max_dist, subclass, -Z ) %>%
    dplyr::select( id, subclass, Y, Z, weights, feasible, adacal, max_dist )

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

                        aa <- estimate_ATT( df_subclass, treatment="Z", outcome="Y" )
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
      make_tx_axis( ggd_att )

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
    make_tx_axis( ggd_att )

  # ylim(c(-0.1, 0.2))

  if ( caliper_plot == "both" ) {
    return( list( plot_SATT = plot_SATT,
                  plot_max_caliper_size = plot_max_caliper_size,
                  table = ggd_att ) )
  }
  plot_SATT

}


