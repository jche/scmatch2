library( CSM )
library(tidyverse)
source("./R/diagnostic_plots.R")
# Test the degree of overlap set by prop_nc_unif
if ( FALSE ){
  sample_dat <-
    gen_one_toy(ctr_dist=0.5,
                prop_nc_unif = 1/3)
  names(sample_dat)
  ggplot( sample_dat, aes( X1, X2, col=as.factor(Z) ) ) +
    geom_point() +
    coord_fixed()
  sample_dat$tau = sample_dat$Y1 - sample_dat$Y0
  skimr::skim( sample_dat )
}

set.seed(123)
dat_high_overlap <-
  gen_one_toy(ctr_dist = 0.5,
              prop_nc_unif = 2/3)

### Match
mtch <- get_cal_matches( dat_high_overlap,
                         metric = "maximum",
                         scaling = 8,
                         caliper = 1,
                         rad_method = "adaptive",
                         est_method = "scm" )
mtch_table <-
  full_unit_table(mtch, nonzero_weight_only = F )

### Plot for what's happening in the class
plot_subclass_weights <-
  function(mtch_table, subclass) {

  filtered_mtch <- mtch_table %>%
    filter(subclass == !!subclass)

  plot <- filtered_mtch %>%
    mutate(Z = factor(Z, levels = c(FALSE, TRUE), labels = c("Control", "Treatment"))) %>%
    ggplot(aes(x = X1, y = X2, col = weights)) +  # Color by weights
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
    theme_classic() +
    labs(x = TeX("$X_1$"),
         y = TeX("$X_2$"),
         color = "Weights") +
    xlim(c(0,0.5)) +
    ylim(c(0,0.5))

  return(plot)
}


plot_high_overlap <-
  plot_subclass_weights(
    mtch_table, subclass = 4)  # Replace with desired subclass
print(plot_high_overlap)


# ### Next have a plot on weights on all units
# generate_plot <- function(weight_threshold, prop_nc_unif, file_name) {
#
#   set.seed(123)
#
#   # Generate the toy dataset
#   dat_high_overlap <- gen_one_toy(ctr_dist = 0.5, prop_nc_unif = prop_nc_unif)
#
#   # Get matching results
#   mtch <- get_cal_matches(dat_high_overlap,
#                           metric = "maximum",
#                           scaling = 8,
#                           caliper = 1,
#                           rad_method = "adaptive",
#                           est_method = "scm")
#
#   # Create the mtch table and aggregate it
#   mtch_table <- full_unit_table(mtch, nonzero_weight_only = F)
#   mtch_table_agg <- mtch_table %>% agg_co_units()
#
#   # Filter and create the plot
#   plot <- mtch_table_agg %>%
#     filter(weights > weight_threshold, Z == F) %>%
#     mutate(Z = factor(Z, levels = c(FALSE, TRUE), labels = c("Control", "Treatment"))) %>%
#     ggplot(aes(x = X1, y = X2, col = weights)) +  # Color by weights
#     geom_point() +
#     scale_color_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
#     theme_classic() +
#     labs(x = TeX("$X_1$"),
#          y = TeX("$X_2$"),
#          color = "Weights")  # Label the color legend
#
#   # Save the plot to a file
#   ggsave(filename = file_name, plot = plot)
# }
#
# # Loop through different values of prop_nc_unif and save each plot
# for (prop_nc_unif in c(1/3, 2/3, 3/3)) {
#   file_name <- paste0("plot_prop_nc_unif_", prop_nc_unif * 3, ".png")
#   generate_plot(
#     weight_threshold = 0.65,
#     prop_nc_unif = prop_nc_unif,
#     file_name = file_name)
# }
#
#
# #
# #
# #
# # ## Average position for treated
# # mtch_table %>%
# #   group_by(Z) %>%
# #   summarise(across(starts_with("X"), mean, na.rm = TRUE)) %>%
# #   ungroup() %>%
# #   ggplot(aes(x = X1, y = X2, )) +
# #   geom_point()
