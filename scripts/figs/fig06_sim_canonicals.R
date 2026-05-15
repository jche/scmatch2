# Plotting script: scripts/figs/fig_sim_canonicals.R
# full analysis of simulation results
library(tidyverse)
source("scripts/lib/plot_sim.R")

# canonical sims-bias_mse ----------------------------------------------------------
res_acic <- read_csv("data/outputs/sims-bias_mse/acic_combined.csv")
res_hain <- read_csv("data/outputs/sims-bias_mse/hainmueller_combined.csv")
res_kang <- read_csv("data/outputs/sims-bias_mse/kang_combined.csv")

res <- bind_rows(res_kang, res_hain, res_acic) %>%
  relocate( sim_type )

table( res$sim_type )


# explore the sim results -----

if ( FALSE ) {

  res_L <- pivot_longer(res_kang,
                        diff:last_col(), names_to="method") %>%
    dplyr::select( -seed, -elapsed_time_secs )
  res_L

  res_L %>% group_by( method ) %>%
    summarise( mean_na = mean( is.na( value ) ),
               SE = sd( value, na.rm=TRUE ),
               sd_tau = sd( true_ATT ) ) %>%
    knitr::kable()

    # Look at the distribution of estimates to see if we have outliers
    ggplot( res_L, aes( method, value ) ) +
      geom_boxplot() +
      coord_flip()

  org_df <- summarize_bias_rmse( res_L )
  org_df %>%
    mutate( del = rmse^2 - se^2 - bias^2 )

  arrange( org_df, rmse )
  ggplot( org_df, aes( bias_mcse, se ) ) +
    geom_point()

  RMSE_plot( res_acic,
            xlab = "Value",
            ylab = "Method", winz = 2
  )

  ggplot( org_df, aes( se, abs(bias) ) ) +
#    geom_point() +
    geom_text( aes( label = method ), size=3 ) +
    theme_minimal()

}


# check rmse --------------------------------------------------------

require(patchwork)

acic_plot_sim_type( "kang",
                    title="Kang & Schafer",
                    xlab="Method",
                    legend.position=c(0.75, 0.25)) +
  acic_plot_sim_type( "hainmueller",
                      ylab="Value", title="Hainmueller") +
  acic_plot_sim_type( "acic",
                      title="ACIC 2016")

ggsave("figures/sim_canonical_results.pdf", width=8, height=3.5, units="in")


# Contour plots: bias vs SE with RMSE contour lines -------------------------

res_toy <- readr::read_csv("data/outputs/sims-bias_mse/toy_combined.csv",
                           show_col_types = FALSE)

make_contour_df <- function(res_long) {
  summarize_bias_rmse(res_long) %>%
    mutate(
      group = case_when(
        method == "csm_scm"       ~ "CSM",
    #    method == "causal_forest" ~ "HighD",
    #    method == "twang" ~ "HighD",
    #    method == "or_bart" ~ "Bart",
        method == "cem_avg"       ~ "CEM",
    #    method == "or_lm"        ~ "LM",
    #    str_detect(method, "ipw") ~ "IPW",
        TRUE                      ~ "Other"
      )
    )
}

res_toy_L <- pivot_longer(res_toy, diff:last_col(), names_to = "method")
res_kang_L <- pivot_longer(res_kang, diff:last_col(), names_to = "method")
res_hain_L <- pivot_longer(res_hain, diff:last_col(), names_to = "method")
res_acic_L <- pivot_longer(res_acic, diff:last_col(), names_to = "method")

p_toy <- contour_plot(make_contour_df(res_toy_L),
                       focus = "csm_scm", title = "Toy Example")
p_kang <- contour_plot(make_contour_df(res_kang_L),
                       focus = "csm_scm", title = "Kang & Schafer")
p_hain <- contour_plot(make_contour_df(res_hain_L),
                       focus = "csm_scm", title = "Hainmueller")
p_acic <- contour_plot(make_contour_df(res_acic_L),
                       focus = "csm_scm", title = "ACIC 2016")

p_toy + p_kang + p_hain + p_acic +
  plot_layout(guides = "collect", nrow=2) &
  theme(legend.position = "bottom")

ggsave("figures/sim_canonical_contour.pdf", width = 10, height = 10, units = "in")







# LaTeX performance tables ---------------------------------------------------
# One table per simulation scenario, sorted by RMSE.
# Columns: Method | |Bias| (MCSE) | SE (MCSE) | RMSE (MCSE)

library(xtable)

# Pretty method name lookup
METHOD_LABELS <- c(
  diff          = "Difference in means",
  onenn         = "1-NN matching",
  csm_scm       = "CSM",
  cem_avg       = "CEM (avg)",
  bal1          = "Balancing weights 1",
  bal2          = "Balancing weights 2",
  or_lm         = "OLS (interactions)",
  or_lm_main    = "OLS (main effects)",
  ps_lm         = "PS regression",
  or_bart       = "BART (outcome)",
  ps_bart       = "BART (PS)",
  aipw1         = "AIPW 1",
  tmle1         = "TMLE 1",
  aipw2         = "AIPW 2",
  tmle2         = "TMLE 2",
  causal_forest = "Causal forest",
  twang         = "TWANG",
  kbal          = "KBAL"
)

make_perf_table <- function(res_long, scenario_label, digits = 3,
                            bold_method = "CSM") {
  fmt <- paste0("%.", digits, "f")

  df <- summarize_bias_rmse(res_long) %>%
    arrange(rmse) %>%
    mutate(
      method_label  = METHOD_LABELS[as.character(method)],
      method_label  = ifelse(is.na(method_label), as.character(method), method_label),
      abs_bias      = abs(bias),
      abs_bias_mcse = bias_mcse
    ) %>%
    select(method_label, abs_bias, abs_bias_mcse, se, se_mcse, rmse, rmse_mcse) %>%
    as.data.frame()

  # Pre-format numeric columns as strings so we can inject \textbf{}
  num_cols <- setdiff(names(df), "method_label")
  df[num_cols] <- lapply(df[num_cols], sprintf, fmt = fmt)

  # Bold every cell in the CSM row
  is_bold <- df$method_label == bold_method
  df[is_bold, ]      <- lapply(df[is_bold, , drop = FALSE],
                               function(x) paste0("\\textbf{", x, "}"))

  colnames(df) <- c(
    "Method",
    "|Bias|", "MCSE",
    "SE",     "MCSE ",    # trailing space avoids duplicate col name
    "RMSE",   "MCSE  "
  )

  # All columns are now character — tell xtable not to reformat (digits = 0)
  xt <- xtable(
    df,
    caption = paste0("Simulation performance: ", scenario_label,
                     ". Methods sorted by RMSE. ",
                     "MCSE = Monte Carlo standard error."),
    label   = paste0("tab:sim_", tolower(gsub("[^A-Za-z0-9]", "_", scenario_label))),
    digits  = 0
  )
  align(xt) <- c("l", "l", "r", "r", "r", "r", "r", "r")
  xt
}

tbl_toy  <- make_perf_table(res_toy_L,  "Toy example")
tbl_kang <- make_perf_table(res_kang_L, "Kang \\& Schafer")
tbl_hain <- make_perf_table(res_hain_L, "Hainmueller")
tbl_acic <- make_perf_table(res_acic_L, "ACIC 2016")

tab_dir  <- here::here("tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
tab_path <- file.path(tab_dir, "sim_canonical_performance.tex")

sink(tab_path)
for (tbl in list(tbl_toy, tbl_kang, tbl_hain, tbl_acic)) {
  print(tbl,
        include.rownames       = FALSE,
        booktabs               = TRUE,
        caption.placement      = "top",
        sanitize.text.function = identity)
  cat("\n\n")
}
sink()

cat("Saved LaTeX tables:", tab_path, "\n")


# # check dropped units -----------------------------------------------------
#
# res %>%
#   ggplot() +
#   geom_density(aes(x=ninf)) +
#   geom_density(aes(x=ninf_cem), color="red") +
#   facet_wrap(~sim, scales="free")
#



