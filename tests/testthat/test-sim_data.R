# test-sim_data.R

test_that("gen_toy_covar returns tibble", {
  out <- gen_toy_covar(10, c(0.2,0.8), c(0.2,0.8), SD = 0.1)
  expect_s3_class(out, "tbl_df")
})

test_that("gen_df_adv returns tibble", {
  out <- gen_df_adv(nc = 20, nt = 10)
  expect_s3_class(out, "tbl_df")
})

test_that("gen_toy_covar_k returns tibble", {
  out <- gen_toy_covar_k(10, centers = list(c(0.2,0.2), c(0.8,0.8)), k = 2, sd = 0.1)
  expect_s3_class(out, "tbl_df")
})

test_that("gen_df_adv_k returns tibble", {
  out <- gen_df_adv_k(nc = 20, nt = 10, k = 3)
  expect_s3_class(out, "tbl_df")
})

test_that("gen_one_toy returns tibble", {
  out <- gen_one_toy(k = 2, nc = 20, nt = 10)
  expect_s3_class(out, "tbl_df")
})

test_that("gen_df_hain returns data frame", {
  out <- gen_df_hain(nt = 10, nc = 20)
  head( out )
  expect_s3_class(out, "data.frame")
})

test_that("gen_df_acic returns tibble (skip if package missing)", {
  testthat::skip_if_not_installed("aciccomp2016")
  out <- gen_df_acic(n = 50, p = 5)
  expect_s3_class(out, "tbl_df")
})

test_that("gen_df_kang returns data frame", {
  out <- gen_df_kang(50)
  expect_s3_class(out, "data.frame")
})
