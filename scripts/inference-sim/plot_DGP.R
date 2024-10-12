
library(dplyr)
library(ggplot2)
library(here)
source(here("scripts/inference-sim", "datagen.R"))

set.seed(123)  # For reproducibility
n <- 2000  # Total sample size
beta_c <- 1 # low overlap
sigma <- 1  # Standard deviation of the error term
beta_0 <- 0.3


beta_c_values  <- c(0, 2, 4)
plot_titles <- c("High Overlap", "Mid Overlap", "Low Overlap")
plots <- list()
for (j in seq_along(beta_c_values)){
  beta_c <- beta_c_values[j]

  dat <- generate_dgp(
    n = n,
    beta_c = beta_c,
    beta_0 = beta_0,
    sigma = sigma,
    n_treated_keep = 10
  )

  plot_title <- paste(plot_titles[j],"; Beta_c =", beta_c)

  plots[[j]] <-
    p <- ggplot(
      dat, aes(x = X, y = Y, color = as.factor(Z))) +
    geom_point() +
    labs(title = plot_title, color = "Treatment (W_i)") +
    theme_minimal()

  if (j < 3){
    plots[[j]] <- plots[[j]] + theme(legend.position = "none")
  }
}

combined_plot <- cowplot::plot_grid(
  plots[[1]],plots[[2]],plots[[3]],
  ncol = 3, align = 'h',
  rel_widths = c(1, 1, 1.4)  # Adjust widths if needed
)

print(combined_plot)
ggsave(
  here("scripts/inference-sim/figures",
       "plot_DGP.png"),
       combined_plot,
       height = 5,
       width = 12)

