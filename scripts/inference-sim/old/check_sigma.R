# Load necessary libraries
library(dplyr)
library(ggplot2)
library(CSM)
library(here)
source(here("scripts/inference-sim", "datagen.R"))
source(here("scripts/inference-sim", "dataviz.R"))


set.seed(123)  # For reproducibility
n <- 2000  # Total sample size
beta_c <- 2  # Coefficient for propensity score model
sigma <- 1  # Standard deviation of the error term
beta_0_values <- c(0.3, 1, 2)  # Different beta_0 values to test
caliper_values <- c(0.1, 0.5, 1, 2)


dat <- generate_dgp(n, beta_c, beta_0, sigma)


bias_matrix <- matrix(NA, nrow = length(beta_0_values), ncol = length(caliper_values))
rownames(bias_matrix) <- paste("Beta_0 =", beta_0_values)
colnames(bias_matrix) <- paste("Caliper =", caliper_values)
tictoc::tic()
for (i in seq_along(beta_0_values)) {
  for (j in seq_along(caliper_values)) {
    # i <- 1; j <- 1
    beta_0 <- beta_0_values[i]
    caliper <- caliper_values[j]

    dat <- generate_dgp(n, beta_c, beta_0, sigma)

    # Matching
    mtch <- get_cal_matches(
      dat,
      covs = "X",
      treatment = "Z",
      metric = "maximum",
      scaling = 1,
      caliper = caliper,
      rad_method = "adaptive",
      est_method = "scm"
    )

    res <- CSM:::get_se_AE(
      mtch,
      outcome = "Y",
      treatment = "Z",
      var_weight_type = "uniform"
    )

    bias <- res$sigma_hat - sigma

    # Store the bias in the matrix
    bias_matrix[i, j] <- bias
  }
}
tictoc::toc()

print(signif(bias_matrix,2))

saveRDS(bias_matrix,
        here::here("data/outputs/inference-sim", "bias_matrix.rds"))
