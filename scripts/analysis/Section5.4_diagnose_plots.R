library(tidyverse)
require(mvtnorm)

library( CSM )
source("./R/diagnostic_plots.R")

produce_histogram_w_j <- function(prop_nc_unif, xlim_range, ylim_range) {
  set.seed(123)
  df_dgp <-
    gen_one_toy(ctr_dist=toy_ctr_dist,
                prop_nc_unif=prop_nc_unif) %>%
    mutate(Y0_denoised = Y0 - noise,
           Y1_denoised = Y1 - noise,
           Y_denoised = Y - noise)
  scaling <- 8

  mtch <- get_cal_matches(
    df_dgp,
    metric = "maximum",
    scaling = scaling,
    caliper = 1,
    rad_method = "adaptive",
    est_method = "scm"
  )

  full_unit_table_mtch <-
    full_unit_table(
      mtch,
      nonzero_weight_only = FALSE
    )

  df_matched_controls <-
    full_unit_table_mtch %>%
    filter(Z == 0) %>%
    group_by(id) %>%
    summarise(w_j = sum(weights))

  df_matched_controls_non_zero_weight <- df_matched_controls %>%
    filter(w_j > 0)

  print(paste0("Number of non-zero weight controls: ",
               nrow(df_matched_controls_non_zero_weight))
  )
  print(paste0("Sum of the weights ",
               sum(df_matched_controls_non_zero_weight$w_j))
  )
  print(paste0("Sum of the weights square ",
               signif( sum(df_matched_controls_non_zero_weight$w_j^2),5)
  )
  )

  hist(df_matched_controls_non_zero_weight$w_j,
       main = paste("Histogram of w_j for prop_nc_unif =", prop_nc_unif),
       xlim = xlim_range,
       ylim = ylim_range,
       xlab = "w_j",
       ylab = "Frequency")
}

# Determine the common x and y limits based on all three datasets
xlim_range <- c(0, 10)  # Adjust based on your data
ylim_range <- c(0, 60) # Adjust based on your data

# Generate histograms for prop_nc_unif values 1/3, 2/3, 3/3
par(mfrow = c(1, 3))  # Set up the plotting area for 3 plots side by side
produce_histogram_w_j(1/3, xlim_range, ylim_range)
produce_histogram_w_j(2/3, xlim_range, ylim_range)
produce_histogram_w_j(3/3, xlim_range, ylim_range)
par(mfrow=c(1,1))


produce_histogram_times_matched <- function(prop_nc_unif, xlim_range, ylim_range) {
  # prop_nc_unif <- 1
  set.seed(123)
  df_dgp <-
    gen_one_toy(ctr_dist=toy_ctr_dist,
                prop_nc_unif=prop_nc_unif) %>%
    mutate(Y0_denoised = Y0 - noise,
           Y1_denoised = Y1 - noise,
           Y_denoised = Y - noise)
  scaling <- 8

  mtch <- get_cal_matches(
    df_dgp,
    metric = "maximum",
    scaling = scaling,
    caliper = 1,
    rad_method = "adaptive",
    est_method = "scm"
  )

  full_unit_table_mtch_nonzero_weight <-
    full_unit_table(
      mtch,
      nonzero_weight_only = TRUE
    )


  tmp <-
    full_unit_table_mtch_nonzero_weight %>%
    filter(Z == 0) %>%
    group_by(id) %>%
    summarise(times_of_used = n())

  k <- 1
  ids_times_more_than_k <-
    tmp %>% filter(times_of_used > k)
  create_toy_df_plot(
    dat_low_overlap %>%
      filter(id %in% ids_times_more_than_k$id)
  )

  hist(tmp$times_of_used,
       main = paste("Histogram of times_of_used for prop_nc_unif =", prop_nc_unif),
       xlim = xlim_range,
       ylim = ylim_range,
       xlab = "Times of Used",
       ylab = "Frequency")
}

# Determine the common x and y limits based on all three datasets
xlim_range <- c(0, 20)  # Adjust based on your data
ylim_range <- c(0, 50)  # Adjust based on your data

# Generate histograms for prop_nc_unif values 1/3, 2/3, 3/3
par(mfrow = c(1, 3))  # Set up the plotting area for 3 plots side by side
produce_histogram_times_matched(1/3, xlim_range, ylim_range)
produce_histogram_times_matched(2/3, xlim_range, ylim_range)
produce_histogram_times_matched(3/3, xlim_range, ylim_range)




produce_matched_control_plot <-
  function(prop_nc_unif) {
    # prop_nc_unif <- 1
    set.seed(123)
    df_dgp <-
      gen_one_toy(ctr_dist=toy_ctr_dist,
                  prop_nc_unif=prop_nc_unif) %>%
      mutate(Y0_denoised = Y0 - noise,
             Y1_denoised = Y1 - noise,
             Y_denoised = Y - noise)
    scaling <- 8

    mtch <- get_cal_matches(
      df_dgp,
      metric = "maximum",
      scaling = scaling,
      caliper = 1,
      rad_method = "adaptive",
      est_method = "scm"
    )

    full_unit_table_mtch_nonzero_weight <-
      full_unit_table(
        mtch,
        nonzero_weight_only = TRUE
      )


    tmp <-
      full_unit_table_mtch_nonzero_weight %>%
      filter(Z == 0) %>%
      group_by(id) %>%
      summarise(times_of_used = n())

    k <- 3
    ids_times_more_than_k <-
      tmp %>% filter(times_of_used > k)
    create_toy_df_plot(
      dat_low_overlap %>%
        filter(id %in% ids_times_more_than_k$id)
    )

  }



matched_control_plot_low <-
  produce_matched_control_plot(1/3)
matched_control_plot_mid <-
  produce_matched_control_plot(2/3)
matched_control_plot_high <-
  produce_matched_control_plot(3/3)
matched_control_combined_plot <- cowplot::plot_grid(
  matched_control_plot_low,
  matched_control_plot_mid,
  matched_control_plot_high,
  ncol = 3, align = 'h',
  rel_widths = c(1, 1, 1)  # Adjust widths if needed
)
matched_control_combined_plot
