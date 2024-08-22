# Function to generate DGP
generate_dgp <-
  function(n, beta_c, beta_0, sigma,
           n_treated_keep = 10,
           prop_treated_keep = NULL) {

  # Step 1: Generate X_i ~ N(0, 1)
  X_i <- runif(n, min=-2,max=2)

  # Step 2: Generate W_i ~ Bernoulli(expit(beta_c * X_i))
  expit <- function(z) {
    return(exp(z) / (1 + exp(z)))
  }

  W_i <- rbinom(n, 1, prob = expit(beta_c * X_i))

  # Step 3: Randomly remove p of treated units
  treated_indices <-
    which(W_i == 1)
  if (!is.null(n_treated_keep)){
    n_treated_keep <- n_treated_keep
  }else if (!is.null(prop_treated_keep)){
    n_treated <- length(treated_indices)
    n_treated_keep <- n_treated * prop_treated_keep
  }
  keep_treated_indices <-
    sample(treated_indices,
           size = n_treated_keep)
  final_indices <-
    sort(c(keep_treated_indices, which(W_i == 0)))

  X_i <- X_i[final_indices]
  W_i <- W_i[final_indices]

  # Step 4: Generate Y_j,t based on Y_j,t ~ beta_0 * X_t + epsilon_j,t
  Y_jt <- beta_0 * X_i + rnorm(length(X_i), mean = 0, sd = sigma)

  # Return a data frame with the generated data
  data <- data.frame(X = X_i, Z = W_i, Y = Y_jt)
  return(data)
}
