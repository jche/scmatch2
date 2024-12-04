library(testthat)
library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))

test_that("DGP is proper", {
  set.seed(123)
  dgp_data <- generate_dgp(N1 = 5, N0 = 1000, panel = "A")
})


test_that("match_controls worked well", {
  set.seed(123)
  dgp_data <- generate_dgp(N1 = 5, N0 = 1000, panel = "A")
  matched_pairs <- match_controls(dgp_data$treated, dgp_data$control, M=4)
})

test_that("compute_overlap_statistics worked well",{
  set.seed(123)
  dgp_data <- generate_dgp(N1 = 5, N0 = 1000, panel = "A")
  matched_pairs <- match_controls(dgp_data$treated, dgp_data$control, M=4)

  source(here("scripts/new-inference/utils-replicate-ferman.R"))
  overlap_stats <- compute_overlap_statistics(matched_pairs)
})
