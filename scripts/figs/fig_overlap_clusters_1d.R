# Define the functions
f <- function(x) {
  0.1 * x^2
}

# Linear effect for treated units (c = 0.2)
f_1 <- function(x) {
  f(x) + 0.2 * x
}

# Create a sequence of x values for plotting the line
x_upper <- 4
x_values <- seq(0, x_upper, by = 0.1)

# Compute y values using the functions
y_values <- f(x_values)
y_values_treated <- f_1(x_values)

# Define matched set points
x_points <- c(10, 35, 90, 95, 100) / 100 * x_upper
y_points <- f(x_points)

# Additional modified points
sig <- 0.8
errors <- c(sig, -sig, -sig, sig, -sig)
y_points_modified <- f(x_points) + errors

# Define treated points (one between green points, one between red and purple, one between purple and blue)
x_treated <- c(0.8, 3.7, 3.9)  # Positioned strategically
y_treated_base <- f_1(x_treated)    # Using the treated function

# Add errors to treated points
treated_errors <- c(0.3, -0.3, 0.4)  # Custom errors for each treated point
y_treated <- y_treated_base + treated_errors  # Apply errors

# Create the plot with larger axis labels
par(mar = c(5, 5, 4, 2) + 0.1)  # Increase margins for larger labels

# Plot the control function
plot(x_values, y_values, type = "l", col = "black", lwd = 2,
     xlab = "", ylab = "",
     xlim = c(0, x_upper),
     ylim = c(-1, 3.5),
     main = "Control vs. Treated Units")

# Add axis labels with larger font
title(xlab = "x", ylab = "f(x)", cex.lab = 1.5)

# Add the treated function line
lines(x_values, y_values_treated, col = "red", lwd = 2, lty = 2)

# Add the control points to the plot
points(x_points[1:2], y_points_modified[1:2], col = "green", pch = 19, cex = 1.2)
points(x_points[3], y_points_modified[3], col = "red", pch = 19, cex = 1.2)
points(x_points[4], y_points_modified[4], col = "purple", pch = 19, cex = 1.2)
points(x_points[5], y_points_modified[5], col = "blue", pch = 19, cex = 1.2)

# Add vertical lines to indicate residuals
for (i in seq_along(x_points)) {
  lines(c(x_points[i], x_points[i]),
        c(y_points[i], y_points_modified[i]),
        col = "orange", lwd = 2)
}

# Add treated points as triangles (pch = 17) with different colors
points(x_treated[1], y_treated[1], col = "green", pch = 17, cex = 1.5)  # First treated point in green
points(x_treated[2], y_treated[2], col = "red", pch = 17, cex = 1.5)    # Second treated point in red
points(x_treated[3], y_treated[3], col = "blue", pch = 17, cex = 1.5)   # Third treated point in blue

# Add vertical lines to indicate errors for treated points
for (i in seq_along(x_treated)) {
  lines(c(x_treated[i], x_treated[i]),
        c(y_treated_base[i], y_treated[i]),
        col = "orange", lwd = 2, lty = 2)
}

# Add a simplified legend (only for control and treated units)
legend("topleft",
       legend = c("Control Units", "Treated Units"),
       col = c("black", "black"),
       pch = c(19, 17),
       pt.cex = c(1.2, 1.5),
       cex = 1.2,
       bg = "white",
       box.lty = 1)

ggplot2::ggsave("fig_overlap_clusters_1d.png",
       width = 8,
       height = 6,
       dpi = 300,
       units = "in",
       device = "png",
       path = "figures/")
