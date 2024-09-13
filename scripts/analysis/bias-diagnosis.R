library(CSM)
library(tidyverse)
library(ggplot2)
library(latex2exp)



# Define subclass values and degree of overlap values
subclass_values <- 1:10
prop_nc_unif_values <- setdiff(seq(0, 1, length.out = 10), 0) #exclude 0

# Create directories to store plots and tables if they don't exist
dir.create("figures/bias-diagnosis", recursive = TRUE, showWarnings = FALSE)

# Initialize an empty list to store bias results for all subclasses and overlap values
bias_results <- list()


# Create directories to store the matching tables
dir.create(here::here("data/bias-diagnosis"), recursive = TRUE, showWarnings = FALSE)

# Step 1: Loop over prop_nc_unif and save mtch_table
for (prop_nc_unif in prop_nc_unif_values) {
  # Set seed to make sure treated is the same set of units
  set.seed(123)

  # Generate data for current prop_nc_unif
  dat <- gen_one_toy(ctr_dist = 0.5, prop_nc_unif = prop_nc_unif)

  # Perform matching
  mtch <- get_cal_matches(dat, metric = "maximum", scaling = 8, caliper = 1, rad_method = "adaptive", est_method = "scm")
  mtch_table <- full_unit_table(mtch, nonzero_weight_only = F)

  # Save the matching table for this prop_nc_unif
  saveRDS(mtch_table, paste0("data/bias-diagnosis/mtch_table_overlap_", round(prop_nc_unif, 2), ".rds"))
}

# Initialize an empty list to store bias results
bias_results <- list()


plot_subclass_weights <- function(mtch_table, subclass) {
  # Filter for the specific subclass
  filtered_mtch <- mtch_table %>%
    filter(subclass == !!subclass)

  # Extract the treated unit (Z == 1)
  treated_unit <- filtered_mtch %>%
    filter(Z == 1) %>%
    slice(1)

  # Calculate bias for the current subclass
  bias_subclass <- filtered_mtch %>%
    group_by(Z) %>%
    summarise(Y0_Z = sum(Y0 * weights)) %>%
    summarise(bias = diff(Y0_Z)) %>%
    pull(bias)

  # Extract the position of the treated unit (X1, X2)
  X1_treated <- treated_unit$X1
  X2_treated <- treated_unit$X2

  # Number of matched controls (excluding the treated unit)
  num_controls <- nrow(filtered_mtch) - 1

  # Generate the plot
  plot <- filtered_mtch %>%
    mutate(Z = factor(Z, levels = c(FALSE, TRUE), labels = c("Control", "Treatment"))) %>%
    ggplot(aes(x = X1, y = X2, col = weights)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_classic() +
    labs(
      x = TeX("$X_1$"),
      y = TeX("$X_2$"),
      color = "Weights",
      title = paste0(
        "Treated at (X1, X2): (", round(X1_treated, 2), ", ", round(X2_treated, 2), ") | ",
        "Bias: ", round(bias_subclass, 3), "\n",  # Line break for the second line
        "Matched Controls: ", num_controls
      )
    ) +
    xlim(c(0, 1)) +  # Fixed range for X1
    ylim(c(0, 1))    # Fixed range for X2

  return(plot)
}


# Function to generate and save plots for each subclass
generate_plots <- function(prop_nc_unif_values, subclass_values) {
  for (prop_nc_unif in prop_nc_unif_values) {
    # Load the saved mtch_table for the current prop_nc_unif
    mtch_table <- readRDS(paste0("data/bias-diagnosis/mtch_table_overlap_", round(prop_nc_unif, 2), ".rds"))

    # Loop over subclasses to generate plots
    for (subclass_selected in subclass_values) {
      # Generate the plot for the current subclass
      plot <- plot_subclass_weights(mtch_table, subclass_selected)

      # Create directory to store plot for the current subclass
      subclass_dir <- paste0("figures/bias-diagnosis/subclass_", subclass_selected)
      dir.create(subclass_dir, recursive = TRUE, showWarnings = FALSE)

      # Save the plot as PNG
      plot_path <- paste0(subclass_dir, "/overlap_", round(prop_nc_unif, 2), ".png")
      ggsave(plot_path, plot = plot, width = 6, height = 4)
    }
  }
}

# Call the function to generate plots
generate_plots(prop_nc_unif_values, subclass_values)


# Function to generate bias tables
generate_bias_tables <- function(prop_nc_unif_values, subclass_values) {
  bias_results <- list()

  for (prop_nc_unif in prop_nc_unif_values) {
    # Load the saved mtch_table for the current prop_nc_unif
    mtch_table <- readRDS(paste0("data/bias-diagnosis/mtch_table_overlap_", round(prop_nc_unif, 2), ".rds"))

    # Store subclass biases for calculating overall bias
    subclass_biases_with_noise <- c()
    subclass_biases <- c()

    # Loop over subclasses to calculate biases
    for (subclass_selected in subclass_values) {
      # Calculate bias_with_noise (using Y0)
      bias_with_noise <- mtch_table %>%
        filter(subclass == subclass_selected) %>%
        group_by(Z) %>%
        summarise(Y0_Z = sum(Y0 * weights)) %>%
        summarise(bias_with_noise = diff(Y0_Z)) %>%
        pull(bias_with_noise)

      # Calculate bias (using Y0 - noise)
      bias <- mtch_table %>%
        filter(subclass == subclass_selected) %>%
        group_by(Z) %>%
        summarise(Y0_denoised_Z = sum((Y0 - noise) * weights)) %>%
        summarise(bias = diff(Y0_denoised_Z)) %>%
        pull(bias)

      # Add to the list of subclass biases for overall calculation
      subclass_biases_with_noise <- c(subclass_biases_with_noise, bias_with_noise)
      subclass_biases <- c(subclass_biases, bias)

      # Extract relevant information for the treated unit (Z == 1)
      treated_unit <- mtch_table %>%
        filter(subclass == subclass_selected, Z == 1) %>%
        select(Y0, X1, X2, noise) %>%
        mutate(Y0_denoised = Y0 - noise) %>%
        slice(1)  # Get the treated unit

      # Record results in the bias_results list
      bias_results[[length(bias_results) + 1]] <- tibble(
        subclass_selected = as.character(subclass_selected),
        prop_nc_unif = prop_nc_unif,
        Y0_treated = treated_unit$Y0,
        Y0_denoised_treated = treated_unit$Y0_denoised,
        X1_treated = treated_unit$X1,
        X2_treated = treated_unit$X2,
        bias_with_noise = bias_with_noise,
        bias = bias
      )
    }

    # Calculate overall bias (average bias across all subclasses for this prop_nc_unif)
    overall_bias_with_noise <- mean(subclass_biases_with_noise)
    overall_bias <- mean(subclass_biases)

    # Record overall bias in the bias_results list
    bias_results[[length(bias_results) + 1]] <- tibble(
      subclass_selected = "all",
      prop_nc_unif = prop_nc_unif,
      Y0_treated = NA,
      Y0_denoised_treated = NA,
      X1_treated = NA,
      X2_treated = NA,
      bias_with_noise = overall_bias_with_noise,
      bias = overall_bias
    )
  }

  # Combine all results into a single data frame
  bias_table <- bind_rows(bias_results)

  # Save bias table to CSV
  write_csv(bias_table, "tables/bias_diagnosis_table.csv")

  # Return bias table for further use
  return(bias_table)
}

# Call the function to generate the bias table
bias_table <- generate_bias_tables(prop_nc_unif_values, subclass_values)


# Save bias table to CSV
write_csv(bias_table, "tables/bias_diagnosis_table.csv")

# Print the bias table to check the results
print(bias_table)

# Create directory to store bias plots if it doesn't exist
dir.create("figures/bias-diagnosis/bias-plots", recursive = TRUE, showWarnings = FALSE)

# Function to generate and save plots for each subclass
generate_bias_plots <- function(bias_table) {
  # Get unique subclass_selected values (including "all")
  subclass_values <- unique(bias_table$subclass_selected)

  for (subclass_value in subclass_values) {
    # subclass_value <- 1
    subclass_data <- bias_table %>%
      filter(subclass_selected ==
               subclass_value)

    # Generate the plot for the current subclass
    bias_plot <- ggplot(subclass_data,
                        aes(x = prop_nc_unif,
                            y = abs(bias),
                            color = as.factor(subclass_value))) +
      geom_line() +
      geom_point() +
      labs(title = paste("Bias for Subclass",
                         subclass_value),
           x = "Degree of Overlap (prop_nc_unif)",
           y = "Absolute Bias",
           color = "Subclass") +
      theme_minimal() +
      theme(legend.position = "none")

    # Create the file path to save the plot
    plot_path <- paste0("figures/bias-diagnosis/bias-plots/bias_plot_subclass_", subclass_value, ".png")

    # Save the plot
    ggsave(plot_path,
           plot = bias_plot, width = 8, height = 6)
  }
}

# Call the function to generate the bias plots
generate_bias_plots(bias_table)

