library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))

# Parameters
N0 <- 1000
N1_values <- c(5, 10, 25, 50)
M_values <- c(1, 4, 10)
panels <- c("A", "B", "C", "D", "E")
tau_0 <- 0
alpha <- 0.05
num_replicates <- 1000
max_permutations <- 1000

set.seed(123)

generate_all_dgp_and_matched_table(
    N0 = N0,
    N1_values = N1_values,
    M_values = M_values,
    panels = panels,
    num_replicates = num_replicates,
    verbose = 2)

source(here("scripts/new-inference/utils-replicate-ferman.R"))
results <-
  generate_full_table(
    N0, N1_values, M_values,
    panels, tau_0, alpha,
    num_replicates = 100,
    max_permutations,
    result_path = here("scripts/new-inference/outputs/results_table.rds")
  )


## Develop the testing procedure with our method
# # Read one dataset
# full_matched_table <-
#   read_one_matched_table(
#     N1 = 10,
#     M = 4,
#     i = 1,
#     panel = "D")

# # Get the
# library(CSM)
# treatment = "Z"
# outcome = "Y"
# var_weight_type = "uniform"
#
# weighted_var <- get_pooled_variance(
#   matches_table = full_matched_table,
#   outcome = outcome,
#   treatment = treatment,
#   var_weight_type = var_weight_type
# )
# sigma_hat <- sqrt(weighted_var)

# ## Test of get_plug_in_SE worked
# signs <- sample(c(-1, 1), N_T, replace=T)
# signs <- c(1, 1,1,1,1,1,1,1,-1,1)
# unique_subclasses <- full_matched_table %>%
#   filter(Z == TRUE) %>%
#   distinct(subclass) %>%
#   arrange(subclass) %>% # Ensure a stable ordering
#   mutate(sign = signs)
#
# # Join the signs back to the full table
# full_matched_table_signed <- full_matched_table %>%
#   left_join(unique_subclasses, by = "subclass") %>%
#   mutate(weights = weights * sign) %>%
#   select(-sign)
#
# Ns <- calc_N_T_N_C(full_matched_table)
# Ns_signed <- calc_N_T_N_C(full_matched_table_signed)
#
# Ns$N_C_tilde
# Ns_signed$N_C_tilde
#
# # Step 5: Calculate the plug-in standard error
# SE <- get_plug_in_SE(
#   N_T = Ns_signed$N_T,
#   ESS_C = Ns_signed$N_C_tilde,
#   sigma_hat = sigma_hat
# )


results_path <- here("scripts/new-inference/outputs/results_table.rds")

source(here("scripts/new-inference/utils-replicate-ferman.R"))

# Generate LaTeX for hierarchical rejection rate table
latex_rejection_rate <- create_latex_table_hierarchical(
  results_path,
  "rejection_rate",
  "Rejection Rates by Panel and M",
  "tab:rejection_rates"
)

# Output the LaTeX code
cat(latex_rejection_rate)


# latex_time_used <- create_latex_table_hierarchical(
#   results_path,
#   "time_used",
#   "Time Used by Panel and M",
#   "tab:time_used"
# )

latex_avg_shared_controls <- create_latex_table_hierarchical(
  results_path,
  "avg_shared_controls",
  "Average Shared Controls by Panel and M",
  "tab:avg_shared_controls"
)
cat(latex_avg_shared_controls)

latex_avg_shared_treated <- create_latex_table_hierarchical(
  results_path,
  "avg_shared_treated",
  "Average Shared Treated by Panel and M",
  "tab:avg_shared_treated"
)
cat(latex_avg_shared_treated)
