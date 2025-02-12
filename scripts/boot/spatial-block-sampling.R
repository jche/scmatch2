# Load required libraries
library(geoR)  # For spatial functions
library(fields)  # For visualization

# Parameters for the variogram
theta1 <- 1  # Nugget effect
theta2 <- 2  # Sill
theta3 <- 1  # Range parameter
grid_size <- 30  # Grid size (e.g., 30x30 grid)

# Define the isotropic exponential covariance function
exponential_covariance <- function(h, theta1, theta2, theta3) {
  return(theta2 * exp(-theta3 * h) + theta1 * (h == 0))
}

# Generate 2D grid of points
x <- seq(-10, 10, length.out = grid_size)
y <- seq(-15, 15, length.out = grid_size)
grid <- expand.grid(x = x, y = y)

# Compute pairwise distances
dist_matrix <- as.matrix(dist(grid))

# Create covariance matrix using the exponential variogram
cov_matrix <- exponential_covariance(dist_matrix, theta1, theta2, theta3)

# Generate zero-mean stationary Gaussian process
set.seed(123)  # For reproducibility
z <- MASS::mvrnorm(mu = rep(0, nrow(grid)), Sigma = cov_matrix)

# Visualize the process
z_matrix <- matrix(z, nrow = grid_size, ncol = grid_size)
fields::image.plot(x, y, z_matrix, main = "Gaussian Process with Exponential Variogram")
