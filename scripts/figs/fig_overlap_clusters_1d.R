# Define the function f(x)
f <- function(x) {
  0.1 * x^2
}

# Create a sequence of x values for plotting the line
x_upper <- 4
x_values <- seq(0, x_upper, by = 0.1)

# Compute y values using the function
y_values <- f(x_values)

# Define matched set points
x_points <- c(10, 35, 90, 95, 96, 100) /100 * 4  # Updated points to illustrate purple region near 100
y_points <- f(x_points)

# Additional modified points
sig <- 0.8
errors <- c(sig, -sig, -sig, sig, -sig, sig)
y_points_modified <- f(x_points) + errors

# Plot the function
plot(x_values, y_values, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "f(x)",
     xlim = c(0, x_upper),
     ylim = c(-1, 3),
     main = "f(x) = 0.1x * 2")

# Add the points to the plot
points(x_points, y_points, pch = 19)

# Add modified points to the plot
points(x_points[1:2], y_points_modified[1:2], col = "green", pch = 19)
points(x_points[3], y_points_modified[3], col = "red", pch = 19)
points(x_points[4:5], y_points_modified[4:5], col = "purple", pch = 19)
points(x_points[6], y_points_modified[6], col = "blue", pch = 19)
# points(x_points[3:5], y_points_modified[3:5], col = "red", pch = 19)
# points(x_points[4:6], y_points_modified[4:6], col = "blue", pch = 19)

# Add vertical lines to indicate residuals
for (i in seq_along(x_points)) {
  lines(c(x_points[i], x_points[i]),
        c(y_points[i], y_points_modified[i]),
        col = "orange", lwd = 2)
}
