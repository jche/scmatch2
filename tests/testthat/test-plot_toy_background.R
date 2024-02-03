library(testthat)
library(ggplot2)
source("../../analysis/plot_toy.R")
#
# # Create a dummy 'toy_background'
# toy_background <- expand.grid(x1 = 1:5, x2 = 1:5)
# toy_background$z <- runif(nrow(toy_background))
#
# test_that("plot_toy_background returns a ggplot", {
#   p <- plot_toy_background(toy_background)
#   expect_true(is.ggplot(p))
# })
