# Subsampling function
subsampling <- function(data, B, l, statistic_fn) {
  n <- length(data)  # Length of the full dataset
  N <- n - l + 1     # Number of overlapping blocks

  # Generate overlapping blocks
  blocks <- lapply(1:N, function(i) data[i:(i + l - 1)])

  # Compute the statistic on each block
  block_statistics <- sapply(blocks, statistic_fn)

  # Compute the full sample statistic
  full_statistic <- statistic_fn(data)

  # Empirical distribution of the centered and scaled statistic
  scaling_factor <- sqrt(l) / sqrt(n)  # Adjust for block size
  empirical_distribution <- sapply(1:B, function(b) {
    # Randomly sample one block for subsampling
    sampled_block <- sample(1:N, size = 1)
    scaled_centered_stat <- scaling_factor * (block_statistics[sampled_block] - full_statistic)
    return(scaled_centered_stat)
  })

  # Compute subsampling bias and variance
  bias_estimate <- mean(block_statistics) - full_statistic
  variance_estimate <- scaling_factor^2 * var(block_statistics)

  # Return results as a list
  return(list(
    empirical_distribution = empirical_distribution,
    full_statistic = full_statistic,
    bias = bias_estimate,
    variance = variance_estimate
  ))
}

# Example usage with the provided data generation process
set.seed(123)

# Data generation function
generate_data <- function(n) {
  epsilon <- rnorm(n + 1)  # Generate (n + 1) to handle 1-dependence
  X <- epsilon[1:n] + 0.5 * epsilon[2:(n + 1)]
  return(X)
}

# Generate example data
n <- 1000
data <- generate_data(n)

# Define a statistic function (e.g., sample mean)
statistic_fn <- function(data) mean(data)

# Subsampling parameters
B <- 500    # Number of bootstrap samples
l <- 20     # Block size

# Run subsampling
results <- subsampling(data, B, l, statistic_fn)

# Print results
print(results$bias)
print(results$variance)
(true_variance = 3.75 / n)
hist(results$empirical_distribution, main = "Empirical Distribution of Statistic",
     xlab = "Centered and Scaled Statistic", col = "lightblue", border = "white")

## The estimated variance is close to the true variance
##  but it's unclear how close it is.
## Next: make monte-carlo simulation to compare subsampling
##  with other 3 methods in block-bootstrap.R
