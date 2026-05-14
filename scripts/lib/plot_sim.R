
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
prepare_method_comparison_df <- function(df) {
  df %>%
    pivot_longer(diff:last_col(), names_to="method") %>%
    summarize_bias_rmse()
}



my_calc_coverage <- function (data, lower_bound, upper_bound, true_param, criteria = c("coverage",
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
    width <- winsorize(width, winz)
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


my_calc_absolute <- function (data, estimates, true_param, criteria = c("bias", "variance",
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
    estimates <- winsorize(estimates, winz)

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
    dat$stddev_mcse <- sqrt(((K - 1)/K) * sum((sqrt(s_sq_t_j) -
                                                 s_t)^2))
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
summarize_bias_rmse <- function(df_long) {
  prep <- df_long %>%
    filter(method %in% METHODS)

  res <- prep %>%
    group_by(method) %>%
    reframe(
      my_calc_absolute(     estimates  = value,
                            true_param = true_ATT,
                            criteria = c("bias","rmse" ) )
    )

  return( res )

  a = filter( df_long, method == "diff" )
  a
  ggplot( df_long, aes( method, value - true_ATT ) ) +
    geom_boxplot() +
    coord_flip()
  ggplot( a, aes( ))

  bias_res <- prep %>%
    group_by(method) %>%
    reframe(
      simhelpers::calc_bias(pick(everything()),
                            estimates  = value,
                            true_param = true_ATT,
                            na.rm      = TRUE)
    ) %>%
    select(method, bias, bias_mcse)

  rmse_res <- prep %>%
    group_by(method) %>%
    reframe(
      simhelpers::calc_rmse(pick(everything()),
                            estimates  = value,
                            true_param = true_ATT,
                            na.rm      = TRUE)
    ) %>%
    select(method, rmse, rmse_mcse)

  left_join(bias_res, rmse_res, by = "method") %>%
    transmute(
      method,
      Bias      = abs(bias),
      Bias_mcse = bias_mcse,
      RMSE      = rmse,
      RMSE_mcse = rmse_mcse
    ) %>%
    mutate(method = fct_reorder(method, RMSE, min)) %>%
    pivot_longer(
      cols      = c(RMSE, Bias),
      names_to  = "name",
      values_to = "value"
    ) %>%
    mutate(mcse = if_else(name == "RMSE", RMSE_mcse, Bias_mcse)) %>%
    select(method, name, value, mcse)
}



RMSE_plot_outline <- function(org_df,
                              legend.position = "none",
                              show_mcse = TRUE) {
  p <- org_df %>%
    ggplot(aes(y = method, x = value, shape = name))

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



modify_axis_to_RMSE_plot<-
  function(p_outline,title="",
           xlab="",
           ylab=""){
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
                      legend.position="none"){
  org_df <- prepare_method_comparison_df(df)
  p_res <- plot_org_df(org_df,
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
                      legend.position="none") {
  df <- res %>%
    filter(sim == simname)
  return(RMSE_plot(df,title,ylab,xlab,legend.position))
}



acic_plot_sim_type <- function(simname,
                      title="",
                      ylab="",
                      xlab="",
                      legend.position="none") {
  df <- res %>%
    filter(sim_type == simname)
  # if (simname == "hainmueller"){
  #   df <- df %>% select(-twang)
  # }
  return(RMSE_plot(df,title,ylab,xlab,legend.position))
}

