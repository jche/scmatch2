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


## Check the overlap impact
N1 = 10
M = 4
panel = "D"
cat("Running Panel", panel, "N1 =", N1, "M =", M, "\n")
result_has_overlap_adj <-
  compute_rejection_rate(
    N1, N0, M, tau_0, alpha,
    num_replicates=300,
    max_permutations = 200,
    panel)
source(here("scripts/new-inference/utils-replicate-ferman.R"))
result_no_overlap_adj <-
  compute_rejection_rate(
    N1, N0, M, tau_0, alpha,
    num_replicates=300,
    max_permutations=200,
    panel,
    overlap_adj = F)
## Seems no impact -- i think this needs to


## Check the method difference
result_ferman <-
  compute_rejection_rate(
    N1, N0, M, tau_0, alpha,
    num_replicates=300,
    max_permutations = 200,
    panel,
    overlap_adj = T)
source(here("scripts/new-inference/utils-replicate-ferman.R"))
result_pooled <-
  compute_rejection_rate(
    N1, N0, M, tau_0, alpha,
    num_replicates=300,
    max_permutations = 200,
    panel,
    overlap_adj = T,
    sd_method = "pooled")
# This runs 1 hour
saveRDS(result_ferman, here("scripts/new-inference/outputs/result_ferman.rds"))
saveRDS(result_pooled, here("scripts/new-inference/outputs/result_pooled.rds"))

result_pooled <- readRDS(here::here("scripts/new-inference/outputs/result_pooled.rds"))
result_ferman <- readRDS(here::here("scripts/new-inference/outputs/result_ferman.rds"))


### Now i want to run ferman
source(here("scripts/new-inference/utils-replicate-ferman.R"))
N0 = 1000; N1_values = 25; M_values = 10;
panels = "F"; num_replicates = 300;
generate_all_dgp_and_matched_table(
  N0 = 1000,
  N1_values = 25,
  M_values = 10,
  panels = "F",
  num_replicates = 300,
  verbose = 2)
source(here("scripts/new-inference/utils-replicate-ferman.R"))
tau_0 = 0; alpha = 0.1; max_permutations = 100
## How to speed up the procedure? If we know the truth,
#   can we calculate the null distribution only just once?
##  i.e., can we get one critical value?
table_testing_panel_F <-
  generate_full_table(
    N0, N1_values, M_values,
    panels, tau_0, alpha,
    num_replicates = 100,
    max_permutations = 50,
    result_path = here("scripts/new-inference/outputs/table_testing_panel_F.rds")
  )
# This is working really well -- 0.11



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
