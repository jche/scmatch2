
# Script to generate figure of local mathing and residuals, for
# motivating inference strategy

# Define the function f(x)
f <- function(x) {
  0.1 * x^2
}

# Create a sequence of x values for plotting the line
x_upper <- 4
x_values <- seq(0, x_upper, by = 0.1)

# Compute y values using the function
y_values <- f(x_values)

# Initial points to plot
x_points <- c(1, 1.5, 3, 3.2, 3.3)
y_points <- f(x_points)

# Additional modified points
sig=0.5
errors <- c(sig, -sig, sig, -sig, sig)
y_points_modified <-
  f(x_points) + errors

# Plot the line
plot(x_values, y_values, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "Y",
     xlim = c(0.8, x_upper),
     ylim=c(-1,3),
     main = "f(x) = 0.1x^2")

# Add the initial points to the plot
points(x_points, y_points)


# Add the modified points to the plot
points(x_points[1:2],
       y_points_modified[1:2],
       col = "green",
       pch = 19)

points(x_points[3:5],
       y_points_modified[3:5],
       col = "purple",
       pch = 19)

# Add treatment points:


for (i in 1:length(x_points)) {
  lines(c(x_points[i], x_points[i]), c(y_points[i], y_points_modified[i]), col = "orange", lwd = 2)
}

