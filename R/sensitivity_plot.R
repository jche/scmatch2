
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
sensitivity_plot <- function( res,
                              data,
                              min_cal = 0.05,
                              max_cal = 5,
                              R = 30,
                              return_table = FALSE ) {

  cals = seq( min_cal, max_cal, length.out = R )

  args = attr( res, "settings" )
  covs = attr( res, "covariates" )

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

  atts <- map_dfr( sen, get_ATT_estimate ) %>%
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

