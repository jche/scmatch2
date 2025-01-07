set.seed(123)

# Example data generation
# X_n = epsilon_n + 0.5 * epsilon_{n+1}, where epsilon_n ~ N(0, 1)
generate_data <- function(n) {
  epsilon <- rnorm(n + 1)  # Generate (n + 1) to handle 1-dependence
  X <- epsilon[1:n] + 0.5 * epsilon[2:(n + 1)]
  return(X)
}

# Naive Bootstrap
naive_bootstrap <- function(data, B, verbose = FALSE) {
  n <- length(data)
  bootstrap_means <- replicate(B, {
    # if (verbose) cat("Naive Bootstrap: Resample\n")
    sample_data <- sample(data, size = n, replace = TRUE)
    mean(sample_data)
  })
  return(bootstrap_means)
}

# Moving Block Bootstrap (MBB)
moving_block_bootstrap <- function(data, B, block_size, verbose = FALSE) {
  n <- length(data)
  num_blocks <- n - block_size + 1
  blocks <- lapply(1:num_blocks, function(i) data[i:(i + block_size - 1)])

  bootstrap_means <- replicate(B, {
    # if (verbose) cat("MBB: Resample blocks\n")
    sampled_blocks <- sample(blocks, size = ceiling(n / block_size), replace = TRUE)
    bootstrap_sample <- unlist(sampled_blocks)[1:n]  # Truncate to original size
    mean(bootstrap_sample)
  })
  return(bootstrap_means)
}

# Block Bootstrap
block_bootstrap <- function(data, B, block_size, verbose = FALSE) {
  n <- length(data)
  num_blocks <- floor(n / block_size)
  blocks <- split(data[1:(num_blocks * block_size)], rep(1:num_blocks, each = block_size))

  bootstrap_means <- replicate(B, {
    # if (verbose) cat("Block Bootstrap: Resample blocks\n")
    sampled_blocks <- sample(blocks, size = num_blocks, replace = TRUE)
    bootstrap_sample <- unlist(sampled_blocks)
    mean(bootstrap_sample)
  })
  return(bootstrap_means)
}

# Monte Carlo Simulation
monte_carlo_simulation <- function(n, B, mc_reps, block_size, verbose = FALSE) {
  true_mean <- 0  # True mean of the data
  naive_coverage <- 0
  mbb_coverage <- 0
  bb_coverage <- 0

  for (i in 1:mc_reps) {
    if (verbose) cat("Monte Carlo iteration:", i, "of", mc_reps, "\n")
    data <- generate_data(n)
    sample_mean <- mean(data)

    # Naive Bootstrap
    naive_means <- naive_bootstrap(data, B, verbose)
    naive_ci <- quantile(naive_means, c(0.025, 0.975))
    naive_coverage <- naive_coverage + (true_mean >= naive_ci[1] && true_mean <= naive_ci[2])

    # Moving Block Bootstrap
    mbb_means <- moving_block_bootstrap(data, B, block_size, verbose)
    mbb_ci <- quantile(mbb_means, c(0.025, 0.975))
    mbb_coverage <- mbb_coverage + (true_mean >= mbb_ci[1] && true_mean <= mbb_ci[2])

    # Block Bootstrap
    bb_means <- block_bootstrap(data, B, block_size, verbose)
    bb_ci <- quantile(bb_means, c(0.025, 0.975))
    bb_coverage <- bb_coverage + (true_mean >= bb_ci[1] && true_mean <= bb_ci[2])
  }

  # Calculate coverage ratios
  naive_coverage <- naive_coverage / mc_reps
  mbb_coverage <- mbb_coverage / mc_reps
  bb_coverage <- bb_coverage / mc_reps

  return(list(naive = naive_coverage, mbb = mbb_coverage, bb = bb_coverage))
}

# Parameters: in future we want to vary n and block_size
n <- 100             # Sample size
B <- 1000            # Number of bootstrap resamples
mc_reps <- 1000      # Number of Monte Carlo repetitions
block_size <- 10     # Block size for MBB and Block Bootstrap

# Run simulation
results <- monte_carlo_simulation(n, B, mc_reps, block_size, verbose = TRUE)

# Save results
output_dir <- here::here("scripts/boot/output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

saveRDS(results, file = file.path(output_dir, "block_boot_coverage_results.rds"))

# Display results
cat("Coverage Ratios (95% CI):\n")
cat("Naive Bootstrap: ", results$naive, "\n")
cat("Moving Block Bootstrap: ", results$mbb, "\n")
cat("Block Bootstrap: ", results$bb, "\n")
