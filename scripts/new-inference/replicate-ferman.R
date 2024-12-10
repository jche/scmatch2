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
    num_replicates,
    max_permutations,
    result_path = here("scripts/new-inference/outputs/results_table.rds")
  )




results_path <- here("scripts/new-inference/outputs/results_table-8-Dec.rds")

source(here("scripts/new-inference/utils-replicate-ferman.R"))
results_path <- "scripts/new-inference/outputs/results_table.rds"

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
