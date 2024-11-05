
# Generates table that shows how well inference works in the main
# paper.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(tidyr)
library(dplyr)
source( here::here( "scripts/analysis/MCSE_functions.R") )

# load data
R = 500
FNAME =
  here::here(
    paste0("data/outputs/A-E-overlap-by-prop-unif/",
           "A_E_toy_low_mid_high_R=",R,".csv")
  )


res <- read_csv(file =  FNAME)
res

summary( res$N_T )
res %>%
  group_by( deg_overlap ) %>%
  summarise( N_T = mean(N_T),
             N_C_tilde = mean(N_C_tilde) )


res <- res %>%
  mutate( error = att_est - att_true,
          noise = att_est - att_true - bias,
          covered = CI_lower <= att_true & att_true <= CI_upper,
          covered_with_true_SE =
            CI_lower_with_true_SE <= att_true & att_true <= CI_upper_with_true_SE )

table( res$deg_overlap )


res <- res %>%
  mutate(
    bias_and_var = att_est - att_true,
    deg_overlap = factor(
      deg_overlap,
      # levels = c("low", "mid", "high")
      levels =  c("very_low", "low", "mid", "high", "very_high")
    )
  )


table_section_5_4 <- res %>%
  group_by(deg_overlap) %>%
  summarise(
    SE_est = mean(se_AE),
    MCSE_SE_est = MCSE_bias(se_AE),
    # SE_True = mean(true_SE),
    # MCSE_SE_True = MCSE_bias(true_SE),
    SE_True = sd(noise),
    MCSE_SE_True = MCSE_SE(noise),
    # SE_samp_infl = sd(error + bias),
    # MCSE_SE_samp_infl = MCSE_SE(error + bias),
    SE_pop = sd(att_est),
    MCSE_SE_pop = MCSE_SE(att_est),
    N_C_tilde = mean( N_C_tilde ),
    mean_bias = mean(bias),
    MCSE_mean_bias = MCSE_bias(bias),
    coverage = mean(covered),
    MCSE_coverage = MCSE_bias(covered),
    coverage_with_true_SE = mean(covered_with_true_SE),
    RMSE = sqrt( mean( (att_est - att_true)^2 ) ),
    MCSE_RMSE = MCSE_SE( att_est - att_true )
  ) %>%
  arrange(deg_overlap)

table_section_5_4

# make tables directory if it doesn't exist
if( !dir.exists( here::here( "tables" ) ) ){
  dir.create( here::here( "tables" ) )
}

write.csv(table_section_5_4,
          file = here::here( "tables/table_section_5_4.csv") )

# pull all columns with MCSE in their name
mcse = table_section_5_4 %>%
  select(contains("MCSE"))

mcse_summary <- mcse %>%
  summarise_all(list(mean = mean, sd = sd)) %>%
  pivot_longer(
    everything(),
    names_to = c("measure", "stat"),
    names_pattern = "^MCSE_(.*)_(mean|sd)$",
    values_to = "value"
  ) %>%
  pivot_wider( names_from = stat, values_from = value )
mcse_summary

cat( "\n\nMCSE summary:\n" )
print( knitr::kable( mcse_summary, digits=3 ) )

# drop those from the table
table_section_5_4 <- table_section_5_4 %>%
  select(-all_of(colnames(mcse)))


table_section_5_4 <- table_section_5_4 %>%
  dplyr::select( deg_overlap, N_C_tilde, SE_est, SE_True, mean_bias, RMSE, coverage,
                 coverage_with_true_SE )

table_section_5_4$deg_overlap = c( "Very Low", "Low", "Medium", "High", "Very High" )

# KNITR TABLE
cat("\n\nknitr version of table for paper:\n" )

print( knitr::kable( table_section_5_4, digits=2 ) )
names(table_section_5_4 )


xt <- table_section_5_4

colnames( xt ) <- c(  "Overlap", "$\\tilde{N}_C$",
                      "avg \\widehat{SE}", "SE",
                      "Mean(Bias)", "RMSE", "Coverage" )

# Print a latex version of the table
library( xtable )
xt <- xtable(
  xt,
  caption = "Performance of the standard error estimates for the average treatment effect (ATE) and the true standard error (SE) across different degrees of overlap.",
  label = "tab:SE_estimates",
  digits = c(0,0,1,2,2,2,2,2) )

cat("\n\nlatex version of table for paper:\n" )

print( xt, include.rownames = FALSE, include.colnames = TRUE,
       caption.placement = "bottom", sanitize.colnames.function = identity )

