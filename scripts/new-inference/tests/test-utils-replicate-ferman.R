library(testthat)
library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))

test_that("DGP is proper", {
  set.seed(123)
  dgp_data <- generate_dgp(N1 = 5, N0 = 1000, panel = "A")
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
    covs = "X",
    treatment = "Z",
    scaling = 1,
    metric = "euclidean",
    rad_method = "knn",
    k = k
  )

  full_matched_table <- full_unit_table(mtch, nonzero_weight_only = FALSE)
  matched_matrix <- get_matched_control_ids(full_matched_table)

  expect_equal(dim(matched_matrix), c(N1, k),
               info = sprintf("Expected matched_matrix to have dimensions (%d, %d)", N1, k))
})
