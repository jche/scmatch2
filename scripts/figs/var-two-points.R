# Define the function f(x)
f <- function(x) {
  x^2
}

# Create a sequence of x values for plotting the line
x_upper <- 100
x_values <- seq(0, x_upper, by = 0.1)

# Compute y values using the function
y_values <- f(x_values)

# Define matched set points
x_points <- c(10, 35, 95, 96, 100)  # Updated points to illustrate purple region near 100
y_points <- f(x_points)

# Additional modified points
sig <- 5
errors <- c(sig, -sig, sig, -sig, sig)
y_points_modified <- f(x_points) + errors

# Plot the function
plot(x_values, y_values, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "f(x)",
     xlim = c(0, x_upper+5),
     ylim = c(-2000, f(x_upper) + 50),
     main = "Comparison of Lipschitz and Derivative Control")

# Add the points to the plot
points(x_points, y_points, pch = 19)

# Add modified points to the plot
points(x_points[1:2], y_points_modified[1:2], col = "green", pch = 19)
points(x_points[3:5], y_points_modified[3:5], col = "purple", pch = 19)

# Add vertical lines to indicate residuals
for (i in seq_along(x_points)) {
  lines(c(x_points[i], x_points[i]),
        c(y_points[i], y_points_modified[i]),
        col = "orange", lwd = 2)
}

# Add slope annotations
x_green_mid   <- (10 + 35)/2
x_purple_mid  <- (95 + 100)/2
slope_green   <- 2 * x_green_mid
slope_purple  <- 2 * x_purple_mid
arrows(x_green_mid, f(x_green_mid), x_green_mid, f(x_green_mid) + 1000,
       length = 0.1, col = "darkgreen", lwd = 2)
text(x_green_mid, f(x_green_mid) + 1000, labels = paste("Slope ~", round(slope_green,2)),, col = "darkgreen")

arrows(x_purple_mid, f(x_purple_mid), x_purple_mid, f(x_purple_mid) -1000,
       length = 0.1, col = "purple", lwd = 2)
text(x_purple_mid, f(x_purple_mid) - 1000, labels = paste("Slope ~", round(slope_purple,2)),, col = "purple")

# Highlight matched set radii
# Green region
text(x_green_mid, -1200, expression(paste(Delta, "x = 15")), col = "darkgreen")
arrows(10, -500, 35, -500, length = 0.1, code = 3, col = "darkgreen")

# Purple region
text(x_purple_mid, -1200, expression(paste(Delta, "x = 5")), col = "purple")
arrows(95, -500, 100, -500, length = 0.1, code = 3, col = "purple")

