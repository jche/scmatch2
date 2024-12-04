library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))

# Parameters
N0 <- 1000
N1_values <- c(5, 10, 25, 50)
M_values <- c(1, 4, 10)
panels <- c("A", "B", "C", "D", "E")
tau_0 <- 0
alpha <- 0.05
num_replicates <- 2
max_permutations <- 1000

# Run full table generation
set.seed(123)
results <-
  generate_full_table(
    N0, N1_values, M_values,
    panels, tau_0, alpha,
    num_replicates,
    max_permutations,
    result_path = here("scripts/new-inference/outputs/results_table.rds"))

cat("Results saved to 'results_table.rds'.\n")


results_path <- here("scripts/new-inference/outputs/results_table.rds")

# Generate LaTeX for different tables
latex_rejection_rate <- create_latex_table(
  results_path,
  "rejection_rate",
  "Rejection Rates by Panel and M",
  "tab:rejection_rates"
)
cat(latex_rejection_rate)
cat("\n\n")
