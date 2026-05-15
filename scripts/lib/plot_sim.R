
# Analyze results for the main simulations on performance of all the
# various matching approaches.

METHODS <- c("diff", "onenn", "csm_scm", "cem_avg", "bal1", "bal2",
             "or_lm", "or_lm_main", "ps_lm",
             "or_bart", "ps_bart",
             "aipw1", "tmle1",
             "aipw2", "tmle2",
             "causal_forest", "twang", "kbal")



#' Prepare Method Comparison Data Frame
#'
#' Transforms the input data frame, which includes estimates from 14 different
#' methods for generated DGP runs with true ATT, into a long format. Each row
#' in the resulting data frame represents a method's estimate for a run. This
#' function is intended as a preliminary step before calculating bias and RMSE.
#'
#' @param df A data frame containing the runs for a generated DGP with true ATT
#'   and estimates from 14 different methods.
#'
#' @return A long format data frame with each row representing (run x method),
#'   ready for bias and RMSE calculation.
#' @export
#'
#' @examples
#' # Assuming df is your dataset
#' prepared_df <- prepare_method_comparison_df(df)
prepare_method_comparison_df <- function(df, winz = 2) {
  df %>%
    pivot_longer(diff:last_col(), names_to="method") %>%
    summarize_bias_rmse( winz = winz )
}



my_calc_coverage <- function (data, lower_bound, upper_bound, true_param,
                              criteria = c("coverage",
                                           "width"), winz = Inf)
{
  criteria <- match.arg(criteria, choices = c("coverage", "width"),
                        several.ok = TRUE)
  if (!missing(data)) {
    cl <- match.call()
    true_param <- eval(cl$true_param, envir = data, enclos = parent.frame())
    lower_bound <- eval(cl$lower_bound, envir = data, enclos = parent.frame())
    upper_bound <- eval(cl$upper_bound, envir = data, enclos = parent.frame())
  }

  not_miss <- !is.na(lower_bound) & !is.na(upper_bound)
  lower_bound <- lower_bound[not_miss]
  upper_bound <- upper_bound[not_miss]
  K <- length(lower_bound)
  width <- upper_bound - lower_bound
  if (winz < Inf)
    width <- simhelpers:::winsorize(width, winz)
  dat <- tibble::tibble(K_coverage = K)
  if (winz < Inf) {
    dat$width_winsor_pct <- attr(width, "winsor_pct")
    dat$width_winsor_pct_mcse <- sqrt(dat$width_winsor_pct *
                                        (1 - dat$width_winsor_pct)/K)
  }

  if ("coverage" %in% criteria) {
    coverage <- mean(lower_bound <= true_param & true_param <=
                       upper_bound)
    dat$coverage <- coverage
    dat$coverage_mcse = sqrt(coverage * (1 - coverage)/K)
  }
  if ("width" %in% criteria) {
    dat$width <- mean(width)
    dat$width_mcse <- sqrt(var(width)/K)
  }
  return(dat)
}


my_calc_absolute <- function (data, estimates, true_param,
                              criteria = c("bias", "variance",
                                           "stddev", "mse", "rmse"), winz = Inf)
{
  criteria <- match.arg(criteria, choices = c("bias", "variance",
                                              "stddev", "mse", "rmse"), several.ok = TRUE)
  if (!missing(data)) {
    cl <- match.call()
    true_param <- eval(cl$true_param, envir = data, enclos = parent.frame())
    estimates <- eval(cl$estimates, envir = data, enclos = parent.frame())
  }
  true_param <- true_param[!is.na(estimates)]
  estimates <- estimates[!is.na(estimates)]
  if (winz < Inf)
    estimates <- simhelpers:::winsorize(estimates, winz)

  errors = estimates - true_param

  K <- length(errors)
  t_bar <- mean(errors)
  s_t <- sd(errors)
  g_t <- sum((errors - t_bar)^3)/(K * s_t^3)
  k_t <- sum((errors - t_bar)^4)/(K * s_t^4)
  mse <- mean(errors^2)
  t_bar_j <- (K * t_bar - errors)/(K - 1)
  bias_j_sq <- t_bar_j^2
  s_sq_t_j <- ((K - 1) * s_t^2 - (errors - t_bar)^2 * K/(K - 1))/(K - 2)
  rmse_j <- sqrt(bias_j_sq + s_sq_t_j)
  dat <- tibble::tibble(K_absolute = K)
  if (winz < Inf) {
    dat$winsor_pct <- attr(errors, "winsor_pct")
    dat$winsor_pct_mcse <- sqrt(dat$winsor_pct * (1 - dat$winsor_pct)/K)
  }

  if ("bias" %in% criteria) {
    dat$bias <- t_bar
    dat$bias_mcse <- s_t/sqrt(K)
  }
  if ("variance" %in% criteria) {
    dat$var <- s_t^2
    dat$var_mcse <- s_t^2 * sqrt((k_t - 1)/K)
  }
  if ("stddev" %in% criteria) {
    dat$stddev <- s_t
    dat$stddev_mcse <- sqrt(((K - 1)/K) * sum((sqrt(s_sq_t_j) - s_t)^2))
  }
  if ("mse" %in% criteria) {
    dat$mse <- mse
    dat$mse_mcse <- sqrt((1/K) * (s_t^4 * (k_t - 1) + 4 *
                                    s_t^3 * g_t * t_bar + 4 * s_t^2 * t_bar^2))
  }
  if ("rmse" %in% criteria) {
    rmse <- sqrt(mse)
    dat$rmse <- rmse
    dat$rmse_mcse <- sqrt(((K - 1)/K) * sum((rmse_j - rmse)^2))
  }
  return(dat)
}


#' Summarize Bias and RMSE for Each Method
#'
#' Given a long-format data frame where each row represents a method's estimate
#' for a run, this function calculates the bias and RMSE for each method across
#' all runs, together with Monte Carlo standard errors (MCSEs) via the
#' \pkg{simhelpers} package. Methods are reordered by RMSE for easier comparison
#' and visualization.
#'
#' @param df_long A long format data frame with each row representing (run x method),
#'   typically prepared by `prepare_method_comparison_df()`. Must contain columns
#'   \code{method}, \code{value} (estimated ATT), and \code{true_ATT}.
#'
#' @return A long data frame with columns:
#'   \describe{
#'     \item{method}{Factor, reordered by ascending RMSE.}
#'     \item{name}{Character: \code{"RMSE"} or \code{"Bias"}.}
#'     \item{value}{Numeric performance estimate (absolute bias or RMSE).}
#'     \item{mcse}{Monte Carlo standard error for \code{value}. For
#'       \code{"Bias"} this is the MCSE of the signed bias—an approximation
#'       for absolute bias that is exact when bias is clearly non-zero.}
#'   }
#' @export
#'
#' @examples
#' # Assuming df_long is your long-format dataset
#' summary_df <- summarize_bias_rmse(df_long)
#'
summarize_bias_rmse <- function(df_long, winz = 2 ) {
  prep <- df_long %>%
    filter(method %in% METHODS)

  prep %>%
    group_by(method) %>%
    reframe(
      my_calc_absolute(
        estimates  = value,
        true_param = true_ATT,
        criteria   = c("bias", "stddev", "rmse"),
        winz = winz
      )
    ) %>%
    mutate(
      method = fct_reorder(method, rmse, min)
    ) %>%
    rename(
      se = stddev,
      se_mcse = stddev_mcse ) %>%
    select(method, bias, bias_mcse, se, se_mcse, rmse, rmse_mcse)
}



RMSE_plot_outline <- function(org_df,
                              legend.position = "none",
                              show_mcse = TRUE) {
  p <- org_df %>%
    ggplot(aes(y = method, x = value, shape = metric))

  if (show_mcse && "mcse" %in% names(org_df)) {
    p <- p +
      geom_linerange(aes(xmin = value - 2 * mcse,
                         xmax = value + 2 * mcse),
                     linewidth = 0.4, alpha = 0.5)
  }

  p +
    geom_point(size = 1) +
    theme_classic() +
    theme(axis.text.y     = element_text(angle = 0, vjust = 0.5, hjust = 1),
          axis.title.y    = element_text(vjust = -3),
          legend.position = legend.position,
          legend.background = element_rect(linetype = "solid", linewidth = 0.5,
                                           color = "black"))
}



modify_axis_to_RMSE_plot <-  function(p_outline,title="",
                                      xlab="",
                                      ylab="") {
  p_outline +
    scale_y_discrete(labels = c(
      "diff" = "diff",
      "tmle1" = "dr-TMLE1",
      "tmle2" = "dr-TMLE2",
      "aipw1" = "dr-AIPW1",
      "aipw2" = "dr-AIPW2",
      "bal1" = "bal-SBW1",
      "bal2" = "bal-SBW2",
      "or_bart" = "or-BART",
      "or_lm" = "or-LM-int",
      "or_lm_main" = "or-LM",
      "ps_bart" = "ps-BART",
      "ps_lm" = "ps-LM",
      "onenn" = "match-1NN",
      "cem_avg" = "match-CEM",
      "csm_scm" = expression(bold(match-CSM))
    )) +
    scale_shape_manual(values = c(1,2)) +
    labs(title = title,
         y = ylab,
         x = xlab,
         pch = "Metric")
}



RMSE_plot <- function(df,
                      title="",
                      xlab="",
                      ylab="",
                      legend.position="none",
                      table_long = FALSE,
                      winz = 2 ) {

  if ( !table_long ) {
    org_df <- prepare_method_comparison_df(df, winz = winz )
  } else {
    org_df = df
  }

  org_df_long <- org_df %>%
    mutate( bias = abs( bias ) ) %>%
    rename(
      bias_value = bias,
      rmse_value = rmse
    ) %>%
    pivot_longer(
      cols = c(bias_value, bias_mcse, rmse_value, rmse_mcse),
      names_to = c("metric", ".value"),
      names_sep = "_"
    )
  p_res <- plot_org_df(org_df_long,
                       title=title,
                       xlab=xlab,
                       ylab=ylab,
                       legend.position=legend.position)
  return(p_res)

}



plot_org_df <- function(org_df,
                        title="",
                        xlab="",
                        ylab="",
                        legend.position="none"){
  p_outline <-
    RMSE_plot_outline(org_df,
                      legend.position=legend.position)
  p_res <-
    modify_axis_to_RMSE_plot(
      p_outline,
      title,
      xlab,
      ylab)
  return(p_res)
}


acic_plot <- function(simname,
                      title="",
                      ylab="",
                      xlab="",
                      legend.position="none",
                      winz = 2 ) {
  df <- res %>%
    filter(sim == simname)
  return(RMSE_plot(df,title,ylab,xlab,legend.position,
                   winz = winz ))
}



#' Contour plot: SE vs |Bias| with RMSE contour lines
#'
#' Plots each method as a point with SE on the x-axis and |Bias| (or Bias) on
#' the y-axis.  Circular RMSE contour lines (sqrt(bias^2 + se^2) = constant)
#' are drawn in the background.  Points are coloured by method type, and
#' optional ±2 MCSE error bars are drawn on each axis.
#'
#' @param results_df  Tibble with columns: method, bias, bias_mcse, se,
#'   se_mcse, rmse, rmse_mcse, and optionally \code{group}.  Typically the
#'   output of \code{summarize_bias_rmse()}.  If a \code{group} column is
#'   present, points are coloured by it; otherwise all points share one colour.
#' @param abs_bias  If TRUE (default), plot |bias| on the y-axis.
#' @param add_labels  If TRUE (default), add repelled method-name labels.
#' @param add_errorbars  If TRUE (default), draw ±2 MCSE error bars on
#'   both axes.
#' @param focus  Optional character vector of method names.  Each named method
#'   gets a hollow black circle drawn around its point to draw the eye to it.
#'   NULL (default) draws no rings.
#' @param group_colours  Named character vector mapping group labels to colours.
#'   Defaults to \code{CONTOUR_GROUP_COLOURS} which maps CSM→red, CF→blue,
#'   CEM→green, IPW→orange, TMLE→purple, Other→grey70.  Pass your own named
#'   vector to override.
#' @param step  Spacing between RMSE contour lines.  Auto-chosen from the
#'   data range if NULL (default).
#' @param n_contour_steps  Approximate number of contour lines when
#'   \code{step} is auto-chosen (default 5).
#' @param title  Optional plot title string.
#' @return A ggplot object.
#'
# Default colour palette for contour_plot() group levels.
# Override by passing a named character vector to the `group_colours` argument.
CONTOUR_GROUP_COLOURS <- c(
  "CEM"   = "#E41A1C",   # red
  "CF"    = "#377EB8",   # blue
  "CSM"   = "#4DAF4A",   # green
  "IPW"   = "#FF7F00",   # orange
  "LM"   = "#16D6E8",
  "TMLE"  = "#984EA3",   # purple
  "Other" = "grey70",
  "all"   = "grey50"
)

contour_plot <- function(results_df,
                         abs_bias         = TRUE,
                         add_labels       = FALSE,
                         add_errorbars    = TRUE,
                         focus            = NULL,
                         group_colours    = CONTOUR_GROUP_COLOURS,
                         step             = NULL,
                         n_contour_steps  = 5,
                         title            = NULL) {

  df <- results_df

  # If no 'group' column is present, add a dummy so colour mapping still works
  if (!"group" %in% names(df)) {
    df <- df %>% mutate(group = "all")
  }

  # Optionally take absolute value of bias
  if (abs_bias) df <- df %>% mutate(bias = abs(bias))

  # Compute axis limits with a little padding
  max_bias <- max(df$bias, na.rm = TRUE)
  max_se   <- max(df$se,   na.rm = TRUE)
  mm       <- max(max_bias, max_se) * 1.2

  # Auto step: pick a round number that gives ~n_contour_steps lines
  if (is.null(step)) {
    step <- signif(mm / n_contour_steps, 1)
    if (step <= 0) step <- mm / n_contour_steps
  }

  # Background grid for contour lines: z = sqrt(x^2 + y^2) = RMSE
  grid_df <- expand.grid(
    x = seq(0, 1.1 * mm, length.out = 200),
    y = seq(0, 1.1 * mm, length.out = 200)
  )
  grid_df$z <- sqrt(grid_df$x^2 + grid_df$y^2)

  contour_breaks <- seq(step, ceiling(mm / step) * step, by = step)

  p <- ggplot(df, aes(x = se, y = bias)) +
    # RMSE contour lines in background
    geom_contour(
      data        = grid_df,
      aes(x = x, y = y, z = z),
      breaks      = contour_breaks,
      colour      = "#744074",
      alpha       = 0.45,
      linewidth   = 0.4,
      inherit.aes = FALSE
    ) +
    # Points coloured by group
    geom_point(aes(colour = group ), size = 1) +
    coord_fixed(xlim = c(0, mm), ylim = c(0, mm), expand = FALSE) +
    labs(
      x      = "SE",
      y      = if (abs_bias) "|Bias|" else "Bias",
      colour = "Group",
      title  = title
    ) +
    scale_colour_manual(
      values = group_colours,
      # pass-through any group values not in the palette rather than erroring
      na.value = "grey70"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position  = "right"
    )

  # ±2 MCSE error bars on both axes
  if (add_errorbars) {
    p <- p +
      geom_errorbar(
        aes(ymin = bias - 2 * bias_mcse, ymax = bias + 2 * bias_mcse,
            colour = group),
        width = 0, alpha = 0.5, linewidth = 0.5
      ) +
      geom_errorbarh(
        aes(xmin = se - 2 * se_mcse, xmax = se + 2 * se_mcse,
            colour = group),
        height = 0, alpha = 0.5, linewidth = 0.5
      )
  }

  # Black ring around focused methods
  if (!is.null(focus)) {
    focus_df <- df %>% filter(method %in% focus)
    if (nrow(focus_df) > 0) {
      p <- p +
        geom_point(
          data        = focus_df,
          aes(x = se, y = bias),
          shape       = 21,
          colour      = "black",
          fill        = NA,
          size        = 3,
          stroke      = 1,
          inherit.aes = FALSE
        )
    }
  }

  # Repelled method name labels
  if (add_labels) {
    p <- p +
      ggrepel::geom_label_repel(
        aes(label = method, colour = group),
        size          = 2.8,
        box.padding   = 0.5,
        point.padding = 0.2,
        max.overlaps  = Inf,
        show.legend   = FALSE,
        label.size    = 0.2
      )
  }

  p
}


acic_plot_sim_type <- function(simname,
                               title="",
                               ylab="",
                               xlab="",
                               legend.position="none",
                               winz = 2) {
  df <- res %>%
    filter(sim_type == simname)
  # if (simname == "hainmueller"){
  #   df <- df %>% select(-twang)
  # }
  return(RMSE_plot(df,title,ylab,xlab,legend.position,
                   winz = winz ))
}

