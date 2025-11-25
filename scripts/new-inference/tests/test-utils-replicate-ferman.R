library(testthat)
library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))

test_that("DGP is proper", {
  set.seed(123)
  dgp_data <- generate_dgp(N1 = 5, N0 = 1000, panel = "A")
  dgp_data <- generate_dgp(N1 = 25, N0 = 1000, panel = "F")
})


# Note: match_controls is depreciated
# test_that("match_controls worked well", {
#   set.seed(123)
#   dgp_data <- generate_dgp(N1 = 5, N0 = 1000, panel = "A")
#   matched_pairs <- match_controls(dgp_data$treated, dgp_data$control, M=4)
# })

test_that("compute_overlap_statistics worked well",{
  set.seed(123)
  dgp_data <- generate_dgp(N1 = 5, N0 = 1000, panel = "A")
  # matched_pairs <- match_controls(dgp_data$treated, dgp_data$control, M=4)

  # source(here("scripts/new-inference/utils-replicate-ferman.R"))
  # overlap_stats <- compute_overlap_statistics(matched_pairs)
})

test_that("get_matched_control_ids generates correct output", {
  N1 <- 10
  N0 <- 1000
  k <- 4
  panel <- "A"

  dgp_data <- generate_dgp(N1 = N1, N0 = N0, panel = panel)

  mtch <- get_cal_matches(
    df = dgp_data,
    covs = starts_with("X"),
    treatment = "Z",
    scaling = 1,
    metric = "euclidean",
    rad_method = "knn",
    k = k
  )

  full_matched_table <- result_table(mtch, nonzero_weight_only = FALSE)
  matched_matrix <- get_matched_control_ids(full_matched_table)

  expect_equal(dim(matched_matrix), c(N1, k),
               info = sprintf("Expected matched_matrix to have dimensions (%d, %d)", N1, k))
})



test_that("get_matched_control_ids generates correct output", {
  N1 <- 25
  N0 <- 1000
  k <- 10
  panel <- "F"

  dgp_data <- generate_dgp(N1 = N1, N0 = N0, panel = panel)

  mtch <- get_cal_matches(
    df = dgp_data,
    covs = starts_with("X"),
    treatment = "Z",
    scaling = 1,
    metric = "euclidean",
    rad_method = "knn",
    k = k
  )

  full_matched_table <- result_table(mtch, nonzero_weight_only = FALSE)
  matched_matrix <- get_matched_control_ids(full_matched_table)

  expect_equal(dim(matched_matrix), c(N1, k),
               info = sprintf("Expected matched_matrix to have dimensions (%d, %d)", N1, k))
})

## Todo: Develop the testing procedure for the pooled variance
# # Read one dataset
# full_matched_table <-
#   read_one_matched_table(
#     N1 = 10,
#     M = 4,
#     i = 1,
#     panel = "D")

# # Get the
# library(CSM)
# treatment = "Z"
# outcome = "Y"
# var_weight_type = "uniform"
#
# weighted_var <- get_pooled_variance(
#   matches_table = full_matched_table,
#   outcome = outcome,
#   treatment = treatment,
#   var_weight_type = var_weight_type
# )
# sigma_hat <- sqrt(weighted_var)

# ## Test of get_plug_in_SE worked
# signs <- sample(c(-1, 1), N_T, replace=T)
# signs <- c(1, 1,1,1,1,1,1,1,-1,1)
# unique_subclasses <- full_matched_table %>%
#   filter(Z == TRUE) %>%
#   distinct(subclass) %>%
#   arrange(subclass) %>% # Ensure a stable ordering
#   mutate(sign = signs)
#
# # Join the signs back to the full table
# full_matched_table_signed <- full_matched_table %>%
#   left_join(unique_subclasses, by = "subclass") %>%
#   mutate(weights = weights * sign) %>%
#   select(-sign)
#
# Ns <- calc_N_T_N_C(full_matched_table)
# Ns_signed <- calc_N_T_N_C(full_matched_table_signed)
#
# Ns$N_C_tilde
# Ns_signed$N_C_tilde
#
# # Step 5: Calculate the plug-in standard error
# SE <- get_plug_in_SE(
#   N_T = Ns_signed$N_T,
#   ESS_C = Ns_signed$N_C_tilde,
#   sigma_hat = sigma_hat
# )

