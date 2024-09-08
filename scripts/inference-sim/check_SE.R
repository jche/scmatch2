library(dplyr)
library(ggplot2)
library(CSM)
library(here)
library(tictoc)

source(here("scripts/inference-sim", "datagen.R"))
# source(here("scripts/inference-sim", "dataviz.R"))

set.seed(123)  # For reproducibility
n <- 2000  # Total sample size
sigma <- 1  # Standard deviation of the error term
caliper <- 0.1
beta_0 <- 0.3

n_t_values <- c(5, 10, 25, 50)
beta_c_values <- c(0, 2, 4)
R <- 100  # Number of repetitions

# Initialize arrays to store SE estimates and biases
se_estimates <- array(NA, dim = c(length(n_t_values), length(beta_c_values), R, 6))
bias_matrix <- array(NA, dim = c(length(n_t_values), length(beta_c_values), 6))

rownames(bias_matrix) <- paste("N_T =", n_t_values)
colnames(bias_matrix) <- paste("Beta_c (Overlap) =", beta_c_values)

tictoc::tic()

# Simulation over R repetitions
for (r in 1:R) {
  for (i in seq_along(n_t_values)) {
    for (j in seq_along(beta_c_values)) {
      print(paste("r:", r, "i:",i, "j:",j))
      n_t <- n_t_values[i]
      beta_c <- beta_c_values[j]

      dat <- generate_dgp(
        n = n, beta_c = beta_c,
        beta_0 = beta_0,
        sigma = sigma,
        n_treated_keep = n_t,
        prop_treated_keep = NULL
      )

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

      AE_res <- CSM:::get_se_AE(
        mtch,
        outcome = "Y",
        treatment = "Z",
        var_weight_type = "uniform"
      )

      CMSE_true <- CSM:::get_plug_in_SE(AE_res$N_T, AE_res$N_C_tilde, sigma)

      SE_boot_sign <- CSM:::boot_SE(mtch, B = 100, boot_mtd = "sign")
      SE_boot_wild <- CSM:::boot_SE(mtch, B = 100, boot_mtd = "wild")

      # Store SE estimates for each scenario
      se_estimates[i, j, r, 1] <- CMSE_true
      se_estimates[i, j, r, 2] <- AE_res$SE
      se_estimates[i, j, r, 3] <- SE_boot_sign$SE_unif_weight
      se_estimates[i, j, r, 4] <- SE_boot_sign$SE_SCM_weight
      se_estimates[i, j, r, 5] <- SE_boot_wild$SE_unif_weight
      se_estimates[i, j, r, 6] <- SE_boot_wild$SE_SCM_weight
    }
  }
}
tictoc::toc()

saveRDS(
  se_estimates,
  here::here("data/outputs/inference-sim",
             "se_estimates_check_SE.rds"))
