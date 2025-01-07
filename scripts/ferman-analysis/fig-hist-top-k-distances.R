
source( here::here( "scripts/ferman-analysis/02-core-ferman-analysis.R" ) )


proportion_by_UF <- ferman_for_analysis %>%
  group_by(Control, is_sao_paolo) %>%
  summarise(n())


plot_dm <- function(dist_to_plot){
  tibble(d = as.numeric(as.matrix(dist_to_plot))) %>%
    ggplot(aes(d)) +
    geom_histogram(color="black", binwidth=0.03) +
    geom_vline(xintercept = c, col="red") +
    theme_classic() +
    labs(y=NULL,
         x = TeX("$d(X_t, X_j)$")) +
    #xlim(c(0,1)) +
    expand_limits(x=0) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

dist_matrix <-
  data.frame(
    t(ferman_scm$dm_uncapped))
dm_col_sorted <- apply(dist_matrix, 2, sort)

if (FALSE){
  dist_top_3 <- dm_col_sorted[1:3,]
  plot_dm(dist_top_3)
}


# Plot top 1, 2, 3 respectively
library(gridExtra)
dm_top_1 <- dm_col_sorted[1,]
dm_top_2 <- dm_col_sorted[2,]
dm_top_3 <- dm_col_sorted[3,]
plot_dm_all <- gridExtra::grid.arrange(
  plot_dm(dm_top_1),
  plot_dm(dm_top_2),
  plot_dm(dm_top_3)
)

ggsave(
  filename= here::here( "figures/hist-top-k-distances.png" ),
  plot=plot_dm_all,
  width=5.3,
  height=3.7
)
