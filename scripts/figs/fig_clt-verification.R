
# Verify the CLT from the simulated data.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
source( here::here( "scripts/inference-scripts/0_sim_inference_utils.R") )

# Load data
results <- read_csv("data/outputs/4d_bias_inference_Apr2025/combined_results.csv")
results <- results %>%
  mutate(ATT = k * 1.5)

# For each of deg_overlap, compare the histogram of
# sqrt(N_T) * (att_est - bias - att_pop) and compare it with
# Normal(0, V). Show the comparison
# (we need a join of table_section_5_4 and res
# on deg_overlap)

# Merge data with res to match by degree of overlap
results <- results %>%
  mutate(
    standardized_att = (att_est - bias - ATT) / SE
  )


overlap_labels <- c("very_low" = "Very Low Overlap",
                    "low" = "Low Overlap",
                    "mid" = "Medium Overlap",
                    "high" = "High Overlap",
                    "very_high" = "Very High Overlap")

# Plot histograms with Normal(0,1) overlay
ggplot(results, aes(x = standardized_att)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~ deg_overlap, scales = "free", labeller = labeller(deg_overlap = overlap_labels)) +
  labs(
    x = "Standardized ATT Estimate",
    y = "Normalized Frequency (100 MC Replications)"
  ) +
  theme_minimal()

ggsave(here::here("figures", "clt-verification.png"),
       width = 10, height = 6, dpi = 300)

# Result: we successfully get that the distributions are quite similar

# # load data
# R = 500
# FNAME =
#   here::here(
#     paste0("data/outputs/A-E-overlap-by-prop-unif/",
#            "A_E_toy_low_mid_high_R=",R,".csv")
#   )
#
#
# res <- read_csv(file =  FNAME)
# res
#
# summary( res$N_T )
# res %>%
#   group_by( deg_overlap ) %>%
#   summarise( N_T = mean(N_T),
#              N_C_tilde = mean(N_C_tilde) )
#
# res <- res %>%
#   rename( CI_lower = "lower",
#           CI_upper = "upper",
#           CI_lower_with_true_SE = "lower_with_true_SE",
#           CI_upper_with_true_SE = "upper_with_true_SE")
#
# res <- res %>%
#   mutate( error = att_est - att_true,
#           noise = att_est - att_true - bias,
#           covered = CI_lower <= att_true & att_true <= CI_upper,
#           covered_with_true_SE =
#             CI_lower_with_true_SE <= att_true & att_true <= CI_upper_with_true_SE )
#
# table( res$deg_overlap )
#
#
# res <- res %>%
#   mutate(
#     bias_and_var = att_est - att_true,
#     d_pop_r = att_est - att_true,
#     deg_overlap = factor(
#       deg_overlap,
#       levels =  c("very_low", "low", "mid", "high", "very_high")
#     ),
#     true_Var = true_SE^2
#   )
#
#
# table_section_5_4 <- res %>%
#   group_by(deg_overlap) %>%
#   summarise(
#     E_Var = mean(true_Var*N_T) ,
#     att_pop = mean(att_true),
#     V_pop = var(att_true)
#   ) %>%
#   mutate( V = E_Var + V_pop) %>% # TODO: we want one V for each run
#   arrange(deg_overlap)
#
# table_section_5_4
#
#
#
#
# ## Verify through coverage ------
# # Compute the coverage
# # First compute the CI using V
# res <- res %>%
#   mutate(
#     CI_lower_V = (att_est - bias) - 1.96 * sqrt(V/N_T),
#     CI_upper_V = (att_est - bias) + 1.96 * sqrt(V/N_T)
#   ) %>%
#   mutate(
#     covered_V = CI_lower_V <= att_true & att_true <= CI_upper_V,
#   )
# coverage <- res %>% group_by(deg_overlap) %>%
#   summarize(coverage_ratio_V = mean(covered_V))
# coverage
#
#
#
# #### Old thought -----------
# ## This old thought made me realized that we need to
# ##  multiply N_T to Var to get the valid variance
# # # Result: the standardized ATT is too diversely distributed
# # # For example, the sd of standardized_att can be as high as 10
# # sd(res %>% filter(deg_overlap == "low") %>% pull(standardized_att) )
# # mean((res %>% filter(deg_overlap == "low") %>% pull(standardized_att) ))
# # Let's debug this
# # interestingly, if we multiply by * sqrt(V), i.e.,
# # When we compare
# # sqrt(N_T) * (att_est - bias - att_pop) and
# # Normal(0, 1), they are quite similar.
# ggplot(res, aes(x = standardized_att * sqrt(V))) +
#   geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7) +
#   stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
#                 color = "red", linetype = "dashed", size = 1) +
#   facet_wrap(~ deg_overlap, scales = "free") +
#   labs(
#     title = "Comparison of Standardized ATT Estimates with Normal(0,1) Distribution",
#     x = "Standardized ATT Estimate",
#     y = "Density"
#   ) +
#   theme_minimal()
# # in other words, the following sd is quite close to 1
# res %>%
#   mutate(att_T = standardized_att) %>%
#   group_by(deg_overlap) %>%
#   summarise(sd(att_T))
#
# # Result
# # deg_overlap `sd(att_T)`
# # <fct>             <dbl>
# #   1 very_low          1.25
# # 2 low               1.04
# # 3 mid               0.931
# # 4 high              0.928
# # 5 very_high         0.888
#
