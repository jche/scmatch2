METHODS <- c("diff", "onenn", "csm_scm", "cem_avg", "bal1", "bal2",
             "or_lm", "ps_lm",
             "or_bart", "ps_bart",
             "aipw1", "tmle1",
             "aipw2", "tmle2")

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
    pivot_longer(diff:aipw2, names_to="method") %>%
    summarize_bias_rmse()
}

#' Summarize Bias and RMSE for Each Method
#'
#' Given a long-format data frame where each row represents a method's estimate
#' for a run, this function calculates the bias and RMSE for each method across
#' all runs. It reorders methods based on their RMSE values for easier comparison
#' and visualization.
#'
#' @param df_long A long format data frame with each row representing (run x method),
#'   typically prepared by `prepare_method_comparison_df()`.
#'
#' @return A data frame with methods as rows, including their RMSE and Bias,
#'   reordered based on RMSE.
#' @export
#'
#' @examples
#' # Assuming df_long is your long-format dataset
#' summary_df <- summarize_bias_rmse(df_long)
summarize_bias_rmse <- function(df_long) {
  df_long %>%
    filter(method %in% METHODS) %>%
    group_by(method) %>%
    summarize(
      RMSE = sqrt(mean((value-true_ATT)^2)),
      Bias = abs(mean(value-true_ATT)),
    ) %>%
    mutate(method = fct_reorder(method, RMSE, min)) %>%
    pivot_longer(c(RMSE, Bias))
}



RMSE_plot_outline <- function(org_df,
                              legend.position="none"){
  org_df %>%
    ggplot() +
    geom_point(aes(y=method,x=value, shape=name), size=1) +
    theme_classic() +
    theme(axis.text.y = element_text(angle=0, vjust=0.5, hjust=1),
          axis.title.y = element_text(vjust=-3),
          legend.position = legend.position,
          legend.background = element_rect(linetype="solid", linewidth=0.5, color="black"))
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
        "or_lm" = "or-LM",
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
