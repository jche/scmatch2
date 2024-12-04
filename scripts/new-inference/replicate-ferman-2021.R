library(MASS) # For mvrnorm
library(dplyr) # For data manipulation
library(here)

# Function to generate DGP based on the panel
generate_dgp <- function(N1, N0, panel = "A") {
  # Generate covariates
  X_treated <- rnorm(N1, mean = 0, sd = 1)
  X_control <- rnorm(N0, mean = 0, sd = 1)

  # Generate outcomes based on the panel
  if (panel == "A") {
    mu_1 <- function(x) x
    epsilon_1 <- function() rnorm(1, mean = 0, sd = 1)
  } else if (panel == "B") {
    mu_1 <- function(x) (qchisq(pnorm(x), df = 8) - 8) / sqrt(16)
    epsilon_1 <- function() rnorm(1, mean = 0, sd = 1)
  } else if (panel == "C") {
    mu_1 <- function(x) (qchisq(pnorm(x), df = 1) - 1) / sqrt(2)
    epsilon_1 <- function() rnorm(1, mean = 0, sd = 1)
  } else if (panel == "D") {
    mu_1 <- function(x) (qchisq(pnorm(x), df = 1) - 1) / sqrt(2)
    epsilon_1 <- function() (rchisq(1, df = 1) - 1) / sqrt(2)
  } else if (panel == "E") {
    mu_1 <- function(x) (qchisq(pnorm(x), df = 1) - 1) / sqrt(2)
    epsilon_1 <- function() 2 * (rchisq(1, df = 1) - 1) / sqrt(2)
  }

  Y_treated <- sapply(X_treated, function(x) mu_1(x) + epsilon_1())
  Y_control <- rnorm(N0, mean = 0, sd = 1)

  list(
    treated = data.frame(ID = 1:N1, X = X_treated, Y = Y_treated),
    control = data.frame(ID = 1:N0, X = X_control, Y = Y_control)
  )
}

# # Function to implement Ferman's Sign-changes algorithm
# ferman_sign_change_test <- function(dgp_data, M, tau_0, alpha = 0.05, max_permutations = 1000) {
#   treated <- dgp_data$treated
#   control <- dgp_data$control
#
#   # Step 1: Match M nearest neighbors
#   treated <- treated %>%
#     rowwise() %>%
#     mutate(
#       matched_controls = list(
#         control %>%
#           arrange(abs(X - X[ID])) %>%
#           slice_head(n = M) %>%
#           pull(Y)
#       )
#     )
#
#   # Step 2: Compute individual treatment effects
#   treated <- treated %>%
#     mutate(
#       tau_hat = Y - mean(unlist(matched_controls)) - tau_0
#     )
#
#   tau_hat_S <- mean(treated$tau_hat)
#
#   # Step 3: Compute test statistic
#   T_obs <- abs(tau_hat_S) / sd(treated$tau_hat)
#
#   # Step 4: Null distribution
#   S <- treated$tau_hat
#   G <- replicate(max_permutations, sample(c(-1, 1), length(S), replace = TRUE))
#   T_null <- apply(G, 2, function(g) {
#     gS <- g * S
#     mean_gS <- mean(gS)
#     abs(mean_gS) / sqrt(sum((gS - mean_gS)^2) / (length(S) - 1))
#   })
#
#   # Step 5: Decision rule
#   critical_value <- quantile(T_null, probs = 1 - alpha)
#   reject <- as.numeric(T_obs > critical_value)
#
#   reject
# }

match_controls <- function(treated, control, M) {
  matched_pairs <- lapply(1:nrow(treated), function(i) {
    distances <- abs(control$X - treated$X[i])
    matched_indices <- order(distances)[1:M]
    list(
      controls = control$Y[matched_indices],
      control_indices = matched_indices
    )
  })
  return(matched_pairs)
}


ferman_sign_change_test <- function(dgp_data, M, tau_0, alpha = 0.05, max_permutations = 1000) {
  treated <- dgp_data$treated
  control <- dgp_data$control
  N1 <- nrow(treated)

  # Step 1: Match M nearest neighbors and track which controls are matched to which treated
  matched_pairs <- lapply(1:N1, function(i) {
    distances <- abs(control$X - treated$X[i])
    matched_indices <- order(distances)[1:M]
    list(
      controls = control$Y[matched_indices],
      control_indices = matched_indices
    )
  })

  # Calculate individual treatment effects
  tau_i <- sapply(1:N1, function(i) {
    treated$Y[i] - mean(matched_pairs[[i]]$controls) - tau_0
  })

  tau_bar <- mean(tau_i)

  # Test statistic following equation (7) in the paper
  T_obs <- abs(tau_bar) / sqrt(sum((tau_i - tau_bar)^2) / (N1 - 1))

  # Generate sign changes respecting shared neighbors
  shared_neighbors <- matrix(FALSE, N1, N1)
  for(i in 1:N1) {
    for(j in 1:N1) {
      if(i < j) {
        shared <- intersect(matched_pairs[[i]]$control_indices,
                            matched_pairs[[j]]$control_indices)
        shared_neighbors[i,j] <- shared_neighbors[j,i] <- length(shared) > 0
      }
    }
  }

  # Generate valid sign changes
  T_null <- replicate(max_permutations, {
    signs <- sample(c(-1,1), 1)
    for(i in 2:N1) {
      connected <- which(shared_neighbors[i,1:(i-1)])
      if(length(connected) > 0) {
        signs[i] <- signs[connected[1]]
      } else {
        signs[i] <- sample(c(-1,1), 1)
      }
    }
    gS <- signs * tau_i
    gS_bar <- mean(gS)
    abs(gS_bar) / sqrt(sum((gS - gS_bar)^2) / (N1 - 1))
  })

  critical_value <- quantile(T_null, probs = 1 - alpha)
  reject <- as.numeric(T_obs > critical_value)
  reject
}


# Function to perform the experiment for one table entry
compute_rejection_rate <- function(N1, N0, M, tau_0, alpha, num_replicates, max_permutations) {
  rejection_rates <- replicate(num_replicates, {
    dgp_data <- generate_dgp(N1, N0)
    ferman_sign_change_test(dgp_data, M, tau_0, alpha, max_permutations)
  })
  mean(rejection_rates)
}

# Function to compute rejection rates for all combinations
generate_full_table <- function(N0, N1_values, M_values, panels, tau_0, alpha, num_replicates, max_permutations) {
  results <- list()
  time_results <- list()

  for (panel in panels) {
    for (N1 in N1_values) {
      for (M in M_values) {
        cat("Running Panel", panel, "N1 =", N1, "M =", M, "\n")

        # Timing start
        tstart <- proc.time()

        # Compute rejection rate
        rejection_rate <- replicate(num_replicates, {
          dgp_data <- generate_dgp(N1, N0, panel)
          ferman_sign_change_test(dgp_data, M, tau_0, alpha, max_permutations)
        })

        # Timing end
        tend <- proc.time()

        # Store results
        results[[paste(panel, N1, M, sep = "_")]] <- mean(rejection_rate)
        time_results[[paste(panel, N1, M, sep = "_")]] <- tend["elapsed"] - tstart["elapsed"]
      }
    }
  }

  # Save results
  saveRDS(results, here("scripts/new-inference/outputs/rejection_rates.rds") )
  saveRDS(time_results, here("scripts/new-inference/outputs/time_results.rds") )

  list(rejection_rates = results, time_results = time_results)
}

# # Parameters
# N1 <- 5
# N0 <- 1000
# M <- 1
# tau_0 <- 0
# alpha <- 0.10
# num_replicates <- 5000
# max_permutations <- min(1000, 2^N1)
#
# # Compute rejection rate
# {
#   set.seed(123)
#   tstart <- proc.time()
#   rejection_rate <-
#     compute_rejection_rate(N1, N0, M, tau_0, alpha, num_replicates, max_permutations)
#   tend <- proc.time()
# }
#
# # Output the result
# cat("Average Rejection Rate for N1 =", N1, ":", rejection_rate, "\n")
# cat("Time used for running", num_replicates, "replications is", tend["elapsed"] - tstart["elapsed"], "seconds\n")
#

# Parameters
N0 <- 1000
N1_values <- c(5, 10, 25, 50)
M_values <- c(1, 4, 10)
panels <- c("A", "B", "C", "D", "E")
tau_0 <- 0
alpha <- 0.10
num_replicates <- 2
max_permutations <- 1000

# Generate the full table
set.seed(123)
results <- generate_full_table(N0, N1_values, M_values, panels, tau_0, alpha, num_replicates, max_permutations)

# Output results
cat("Rejection rates and time results saved as .rds files.\n")
