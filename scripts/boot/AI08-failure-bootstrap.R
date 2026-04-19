# Bootstrap Failure for Matching Estimators
# Based on Abadie & Imbens (2008)

set.seed(42)

# Load required libraries
library(ggplot2)
library(gridExtra)

# Simulation parameters
n_mc <- 1000      # Number of Monte Carlo simulations
n_boot <- 500     # Number of bootstrap samples
N1 <- 40          # Number of treated units
N0 <- 160         # Number of control units
alpha <- N1 / N0  # Ratio of treated to control
tau_true <- 1.0   # True treatment effect

# Function to generate data under the DGP from Abadie & Imbens (2008)
generate_data <- function(N1, N0, tau) {
  N <- N1 + N0

  # Generate covariates uniformly on [0, 1]
  X <- runif(N)

  # Treatment assignment
  W <- c(rep(1, N1), rep(0, N0))

  # Generate outcomes
  Y <- numeric(N)
  Y[W == 1] <- tau                    # Treated units have outcome = tau
  Y[W == 0] <- rnorm(N0, mean = 0, sd = 1)  # Control units ~ N(0, 1)

  return(list(X = X, W = W, Y = Y))
}

# Function to perform nearest neighbor matching with replacement
nearest_neighbor_matching <- function(X, W, Y) {
  treated_idx <- which(W == 1)
  control_idx <- which(W == 0)

  X_treated <- X[treated_idx]
  X_control <- X[control_idx]

  # Find nearest control for each treated
  matches <- sapply(X_treated, function(x) {
    distances <- abs(x - X_control)
    which.min(distances)
  })

  # Calculate treatment effect
  Y_treated <- Y[treated_idx]
  Y_matched <- Y[control_idx[matches]]

  tau_hat <- mean(Y_treated - Y_matched)

  # Count how many times each control is used as a match
  K <- tabulate(matches, nbins = length(control_idx))

  return(list(tau_hat = tau_hat, K = K, matches = matches))
}

# Function to perform bootstrap for matching estimator
bootstrap_matching <- function(X, W, Y, n_boot) {
  N1 <- sum(W)
  N0 <- sum(1 - W)

  treated_idx <- which(W == 1)
  control_idx <- which(W == 0)

  tau_boot <- numeric(n_boot)

  for (b in 1:n_boot) {
    # Resample treated and control separately
    boot_treated <- sample(treated_idx, size = N1, replace = TRUE)
    boot_control <- sample(control_idx, size = N0, replace = TRUE)
    boot_idx <- c(boot_treated, boot_control)

    X_boot <- X[boot_idx]
    W_boot <- W[boot_idx]
    Y_boot <- Y[boot_idx]

    result <- nearest_neighbor_matching(X_boot, W_boot, Y_boot)
    tau_boot[b] <- result$tau_hat
  }

  return(tau_boot)
}

# Function to calculate analytic variance estimator from Abadie & Imbens (2006)
analytic_variance <- function(X, W, Y, tau_hat) {
  N1 <- sum(W)
  N0 <- sum(1 - W)

  treated_idx <- which(W == 1)
  control_idx <- which(W == 0)

  X_treated <- X[treated_idx]
  X_control <- X[control_idx]

  # Find matches
  matches <- sapply(X_treated, function(x) {
    distances <- abs(x - X_control)
    which.min(distances)
  })

  # Count K_i (number of times control i is used as match)
  K <- tabulate(matches, nbins = N0)

  # For K_sq, count with weights
  K_sq <- numeric(N0)
  for (i in 1:length(matches)) {
    K_sq[matches[i]] <- K_sq[matches[i]] + 1
  }

  # Estimate conditional variance using within-group matching
  sigma_sq <- numeric(length(X))

  # For treated units
  for (i in treated_idx) {
    dist_to_treated <- abs(X[i] - X[treated_idx])
    dist_to_treated[treated_idx == i] <- Inf
    closest_treated <- treated_idx[which.min(dist_to_treated)]
    sigma_sq[i] <- 0.5 * (Y[i] - Y[closest_treated])^2
  }

  # For control units
  for (i in control_idx) {
    dist_to_control <- abs(X[i] - X[control_idx])
    dist_to_control[control_idx == i] <- Inf
    closest_control <- control_idx[which.min(dist_to_control)]
    sigma_sq[i] <- 0.5 * (Y[i] - Y[closest_control])^2
  }

  # Calculate variance components
  Y_treated <- Y[treated_idx]
  Y_matched <- Y[control_idx[matches]]

  var1 <- mean((Y_treated - Y_matched - tau_hat)^2)
  var2 <- sum((K^2 - K_sq) * sigma_sq[control_idx]) / N1^2

  return(var1 / N1 + var2)
}

# Run Monte Carlo simulation
cat("Running Monte Carlo simulation...\n")
tau_estimates <- numeric(n_mc)
boot_vars <- numeric(n_mc)
analytic_vars <- numeric(n_mc)

for (mc in 1:n_mc) {
  if (mc %% 100 == 0) {
    cat(sprintf("  MC iteration %d/%d\n", mc, n_mc))
  }

  # Generate data
  data <- generate_data(N1, N0, tau_true)
  X <- data$X
  W <- data$W
  Y <- data$Y

  # Original estimate
  result <- nearest_neighbor_matching(X, W, Y)
  tau_hat <- result$tau_hat
  K <- result$K
  tau_estimates[mc] <- tau_hat

  # Bootstrap variance
  tau_boot <- bootstrap_matching(X, W, Y, n_boot)
  boot_vars[mc] <- var(tau_boot - tau_hat)  # V^{B,I}

  # Analytic variance
  analytic_vars[mc] <- analytic_variance(X, W, Y, tau_hat)
}

# Calculate true variance from Monte Carlo
true_var <- var(tau_estimates)

# Print results
cat("\n")
cat(strrep("=", 60), "\n")
cat("RESULTS\n")
cat(strrep("=", 60), "\n")
cat(sprintf("True tau: %.4f\n", tau_true))
cat(sprintf("Mean estimated tau: %.4f\n", mean(tau_estimates)))
cat("\nVariance Estimates:\n")
cat(sprintf("  True variance (from MC): %.6f\n", true_var))
cat(sprintf("  Mean bootstrap variance: %.6f\n", mean(boot_vars)))
cat(sprintf("  Mean analytic variance:  %.6f\n", mean(analytic_vars)))
cat("\nStandard Errors:\n")
cat(sprintf("  True SE (from MC): %.6f\n", sqrt(true_var)))
cat(sprintf("  Mean bootstrap SE: %.6f\n", sqrt(mean(boot_vars))))
cat(sprintf("  Mean analytic SE:  %.6f\n", sqrt(mean(analytic_vars))))
cat("\nRelative Bias:\n")
cat(sprintf("  Bootstrap: %.2f%%\n", (mean(boot_vars) - true_var) / true_var * 100))
cat(sprintf("  Analytic:  %.2f%%\n", (mean(analytic_vars) - true_var) / true_var * 100))

# Create visualizations
# Plot 1: Distribution of tau estimates
df_tau <- data.frame(tau_estimates = tau_estimates)
p1 <- ggplot(df_tau, aes(x = tau_estimates)) +
  geom_histogram(bins = 50, alpha = 0.7, fill = "lightblue", color = "black") +
  geom_vline(xintercept = tau_true, color = "red", linetype = "dashed",
             size = 1, aes(linetype = "True τ")) +
  geom_vline(xintercept = mean(tau_estimates), color = "blue", linetype = "dashed",
             size = 1, aes(linetype = "Mean estimate")) +
  labs(x = "Treatment Effect Estimate", y = "Frequency",
       title = "Distribution of Matching Estimates") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plot 2: Variance estimates comparison
df_var <- data.frame(
  value = c(boot_vars, analytic_vars),
  method = rep(c("Bootstrap", "Analytic"), each = n_mc)
)
p2 <- ggplot(df_var, aes(x = method, y = value)) +
  geom_boxplot(fill = c("lightcoral", "lightgreen")) +
  geom_hline(yintercept = true_var, color = "red", linetype = "dashed",
             size = 1) +
  annotate("text", x = 1.5, y = true_var * 1.1, label = "True Variance",
           color = "red") +
  labs(x = "", y = "Variance Estimate",
       title = "Variance Estimator Comparison") +
  theme_minimal()

# Plot 3: Scatter plot of bootstrap vs analytic variance
df_scatter <- data.frame(analytic = analytic_vars, bootstrap = boot_vars)
max_val <- max(c(analytic_vars, boot_vars))
p3 <- ggplot(df_scatter, aes(x = analytic, y = bootstrap)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed",
              size = 1) +
  annotate("text", x = max_val * 0.7, y = max_val * 0.8, label = "45° line",
           color = "red") +
  labs(x = "Analytic Variance", y = "Bootstrap Variance",
       title = "Bootstrap vs Analytic Variance") +
  theme_minimal()

# Plot 4: Bias distribution
boot_bias <- boot_vars - true_var
analytic_bias <- analytic_vars - true_var
df_bias <- data.frame(
  bias = c(boot_bias, analytic_bias),
  method = rep(c("Bootstrap", "Analytic"), each = n_mc)
)
p4 <- ggplot(df_bias, aes(x = bias, fill = method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("Bootstrap" = "lightcoral", "Analytic" = "lightgreen")) +
  labs(x = "Bias (Estimated Var - True Var)", y = "Frequency",
       title = "Distribution of Variance Estimation Bias",
       fill = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine all plots
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)

# Save the plot
ggsave("bootstrap_failure_matching.png", combined_plot, width = 12, height = 10, dpi = 300)

cat("\nPlot saved as 'bootstrap_failure_matching.png'\n")

# Additional analysis: Show how bias changes with alpha
cat("\n")
cat(strrep("=", 60), "\n")
cat("THEORETICAL PREDICTIONS\n")
cat(strrep("=", 60), "\n")

# Calculate theoretical limits
limit_var <- 1 + 1.5 * alpha
limit_boot_var <- 1 + 1.5 * alpha * (5 * exp(-1) - 2 * exp(-2)) / (3 * (1 - exp(-1))) + 2 * exp(-1)

cat(sprintf("Alpha (N1/N0): %.4f\n", alpha))
cat(sprintf("Theoretical limit of N1*Var(tau): %.6f\n", limit_var))
cat(sprintf("Theoretical limit of N1*E[V^{B,I}]: %.6f\n", limit_boot_var))
cat(sprintf("Theoretical bias: %.6f\n", limit_boot_var - limit_var))
cat(sprintf("Observed N1*Var(tau): %.6f\n", N1 * true_var))
cat(sprintf("Observed N1*E[V^{B,I}]: %.6f\n", N1 * mean(boot_vars)))
cat(sprintf("Observed bias: %.6f\n", N1 * (mean(boot_vars) - true_var)))
