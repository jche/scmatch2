library(FNN)
library(dplyr)

set.seed(4044440)

# Generate toy data
dat <- tibble::tibble(
  id = 1:503,
  X1 = runif(503),
  X2 = runif(503),
  Z = c(rep(TRUE, 5), rep(FALSE, 503 - 5)), # Treated: first 5 rows
  noise = rnorm(503),
  Y0 = 3 + 2 * X1 + noise,
  Y1 = 3 + 2 * X1 + 4 * X2 + noise,
  Y = ifelse(Z, Y1, Y0)
)

# Separate treated and control groups
treated <- dat %>% filter(Z == TRUE)
control <- dat %>% filter(Z == FALSE)

# Covariates for matching
treated_cov <- treated %>% select(X1, X2)
control_cov <- control %>% select(X1, X2)

# Perform KNN matching
k <- 4
knn_result <- get.knnx(data = control_cov, query = treated_cov, k = k)

# Extract indices of matches and compute matched control outcomes
matched_indices <- knn_result$nn.index
matched_control_outcomes <- control$Y0[matched_indices]

# Create a result data frame
matching_result <- treated %>%
  mutate(
    matched_Y0 = rowMeans(matched_control_outcomes), # Average if k > 1
    ATT = Y - matched_Y0 # Compute ATT for each treated unit
  )

# Output results
print(matching_result)
