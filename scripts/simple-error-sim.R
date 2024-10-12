

# Set seed for reproducibility
set.seed(123)

# Define parameters
N_T <- 100  # Number of time periods
sigma <- 1  # True sigma value
n_t <- rep(3, N_T)  # Number of observations per time period (constant)

# Generate the data
X <- rnorm(N_T, 0, 1)
Y <-
  lapply(1:N_T,
         function(t) X[t] + rnorm(n_t[t], 0, sigma))


# Compute the estimator S^2
e_t <- sapply(1:N_T, function(t) {
  Y_t <- Y[[t]]
  Y_bar_t <- mean(Y_t)
  e_j_t <- Y_t - Y_bar_t
  s2_t <- sum(e_j_t^2) / (n_t[t] - 1)
  return(n_t[t] * s2_t)
})

N_C <- sum(n_t)
S2 <- sum(e_t) / N_C

# Print the results
cat("Estimated sigma^2: ", S2, "\n")
cat("True sigma^2: ", sigma^2, "\n")

# Optional: Plot the distribution of the estimates
library(ggplot2)
sigma2_estimates <- replicate(1000, {
  X <- rnorm(N_T, 0, 1)
  Y <- lapply(1:N_T, function(t) X[t] + rnorm(n_t[t], 0, sigma))
  e_t <- sapply(1:N_T, function(t) {
    Y_t <- Y[[t]]
    Y_bar_t <- mean(Y_t)
    e_j_t <- Y_t - Y_bar_t
    s2_t <- sum(e_j_t^2) / (n_t[t] - 1)
    return(n_t[t] * s2_t)
  })
  N_C <- sum(n_t)
  S2 <- sum(e_t) / N_C
  return(S2)
})

df <- data.frame(estimate = sigma2_estimates)
ggplot(df, aes(x = estimate)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  geom_vline(aes(xintercept = sigma^2), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of sigma^2 Estimates", x = "Estimate", y = "Frequency") +
  theme_minimal()

