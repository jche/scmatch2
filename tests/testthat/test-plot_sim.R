source("../../analysis/plot_sim.R")

res_toy <-
  readr::read_csv("../../sim_toy_results/toy_spaceship7.csv",
                  show_col_types = F)
res_toy <- res_toy %>%
  rename(tmle2="tmle3",
         aipw2="aipw3")

test_that("prepare_method_comparison_df correctly gives the column names ", {
  result <- prepare_method_comparison_df(res_toy)

  expect_equal(colnames(result), c("method","name","value"))
})

test_that("summarize_df_RMSE correctly computes RMSE and Bias", {
  df_long <- data.frame(method = rep(c("diff", "onenn"), each = 2),
                   value = c(1, 2, 1.5, 2.5),
                   true_ATT = c(1.5, 1.5, 2.5, 2.5))
  result <- summarize_df_RMSE(df_long)
  expected_rmse <- sqrt(mean(c((1-1.5)^2, (2-1.5)^2)))
  expect_equal(result %>%
                 filter(name=="RMSE") %>%
                 pull(value) %>%
                 first(), expected_rmse)
})


test_that("Integration test for plotting functions", {
  plot_toy <- RMSE_plot(res_toy,
            title="Toy Example",
            xlab="Value",
            ylab="Method",
            legend.position=c(0.75, 0.25))
  expect_true(is.ggplot(plot_toy))
})


