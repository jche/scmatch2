

source( here::here( "scripts/ferman-analysis/02-core-ferman-analysis.R" ) )

ferman_for_analysis <- ferman_for_analysis %>%
  mutate(`X1 (2007 score)` = y2007,
         `X2 (2008 score)` = y2008,
         `X3 (2009 score)` = y2009,
         `X4 (pct Sao Paolo)` = is_sao_paolo)

covs <- c( "y2007", "y2008", "y2009", "is_sao_paolo" )

covs_names <-c("X1 (2007 score)", "X2 (2008 score)",
         "X3 (2009 score)", "X4 (pct Sao Paolo)")


love_plot(
  res = ferman_scm,
  covs = covs,
  covs_names = covs_names,
  B=NA)

ggsave(
  filename="figures/love-plot-ferman.png",
  width=8.76,
  height=5.33
)
