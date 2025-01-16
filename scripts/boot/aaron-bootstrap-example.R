# ------------------------------------------------------------
# 1) Setup parameters
# ------------------------------------------------------------
n <- 100                # sample size
mu <- rep(0, n)         # mean vector (all zeros for simplicity)
sigma2 <- 1             # diagonal variance
sigma_small2 <- 0.1     # small off-diagonal variance
k <- 2                  # number of off-diagonal bands to fill

# ------------------------------------------------------------
# 2) Construct Toeplitz covariance matrix
# ------------------------------------------------------------
make_toeplitz_sigma <- function(n, sigma2, sigma_small2, k) {
  Sigma <- matrix(0, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      if(i == j) {
        Sigma[i, j] <- sigma2
      } else if(abs(i - j) <= k) {
        Sigma[i, j] <- sigma_small2
      } else {
        Sigma[i, j] <- 0
      }
    }
  }
  return(Sigma)
}

Sigma <- make_toeplitz_sigma(n, sigma2, sigma_small2, k)

# ------------------------------------------------------------
# 3) Generate data from N(mu, Sigma)
# ------------------------------------------------------------
library(MASS)
set.seed(123)
X <- mvrnorm(1, mu, Sigma)  # single realization, X is length n

# In typical situations, you'd do multiple realizations or treat
# X as your observed data from the process.

# ------------------------------------------------------------
# Naive bootstrap for the mean
# ------------------------------------------------------------
B <- 1000  # number of bootstrap replicates
Xbar <- mean(X)

boot_means_naive <- numeric(B)
for(b in 1:B) {
  X_star <- sample(X, size = n, replace = TRUE)
  boot_means_naive[b] <- mean(X_star)
}

# 95% CI via naive bootstrap
ci_naive <- quantile(boot_means_naive, probs = c(0.025, 0.975))

cat("Naive bootstrap CI:", ci_naive, "\n")


# ------------------------------------------------------------
# Parametric bootstrap that respects estimated covariance
# ------------------------------------------------------------
library(MASS)

# Estimate mu and Sigma from data
mu_hat <- mean(X)  # scalar if univariate
Sigma_hat <- var(X)  # naive univariate variance estimate
# For a general multivariate approach, you'd use cov(X_mat) etc.

# But let's pretend we know the Toeplitz structure (or have estimated it).
# We might use Sigma or Sigma_hat (for real data).

B <- 1000
boot_means_param <- numeric(B)

for(b in 1:B) {
  # parametric sample from N(mu_hat, Sigma)
  X_star <- mvrnorm(1, mu = rep(mu_hat, n), Sigma = Sigma)
  boot_means_param[b] <- mean(X_star)
}

# 95% CI via parametric bootstrap
ci_param <- quantile(boot_means_param, probs = c(0.025, 0.975))

cat("Parametric bootstrap CI:", ci_param, "\n")



########
######## MCMC coverage study
########
# ------------------------------------------------------------
# Coverage comparison
# ------------------------------------------------------------
n_sim <- 100  # number of simulation replications
true_mu <- 0   # we set mu=0 above

coverage_naive <- 0
coverage_param <- 0

for(sim_i in 1:n_sim) {
  # Print progress
  cat("Running simulation:", sim_i, "of", n_sim, "\n")

  # Generate data
  X <- mvrnorm(1, mu, Sigma)

  # NAIVE BOOTSTRAP
  Xbar <- mean(X)
  boot_means_naive <- replicate(1000, mean(sample(X, n, replace=TRUE)))
  ci_naive <- quantile(boot_means_naive, probs=c(0.025, 0.975))
  coverage_naive <- coverage_naive +
    as.numeric(ci_naive[1] <= true_mu & true_mu <= ci_naive[2])

  # PARAMETRIC BOOTSTRAP
  # (Pretend we know Sigma exactly or estimate it)
  boot_means_param <- numeric(1000)
  for(b in 1:1000) {
    X_star <- mvrnorm(1, rep(Xbar, n), Sigma)
    boot_means_param[b] <- mean(X_star)
  }
  ci_param <- quantile(boot_means_param, probs=c(0.025, 0.975))
  coverage_param <- coverage_param +
    as.numeric(ci_param[1] <= true_mu & true_mu <= ci_param[2])

  # Save intermediate results for each simulation
  all_results[[sim_i]] <- list(
    sim_id = sim_i,
    X = X,
    Xbar = Xbar,
    ci_naive = ci_naive,
    ci_param = ci_param
  )
}

# Final results
final_results <- list(
  naive_coverage = coverage_naive / n_sim,
  parametric_coverage = coverage_param / n_sim,
  all_simulations = all_results
)

# Save results to file
output_path <- "scripts/boot/output/aaron-bootstrap-example.rds"
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)  # Ensure directory exists
saveRDS(final_results, file = output_path)

cat("Results saved to", output_path, "\n")
