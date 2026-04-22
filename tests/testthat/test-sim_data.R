# test-sim_data.R

test_that("gen_toy_covar returns tibble", {
  out <- gen_toy_covar(10, c(0.2,0.8), c(0.2,0.8), SD = 0.1)
  expect_s3_class(out, "tbl_df")
})

test_that("gen_df_adv returns tibble", {
  out <- gen_df_adv(nc = 20, nt = 10)
  expect_s3_class(out, "tbl_df")
  expect_equal( nrow(out), 30 )
  expect_equal( sum(out$Z), 10 )

  out <- gen_df_adv(nc = 20, nt = 10, )
  expect_s3_class(out, "tbl_df")
  expect_equal( nrow(out), 30 )
  expect_equal( sum(out$Z), 10 )
  expect_equal( out$Y0, out$Y1 - 1 )

  # Variance function for noise
  f0_fun <- function( x1, x2 ) {
    -100 + x1 * 100
  }
  out <- gen_df_adv(nc = 20, nt = 10, f0_fun = f0_fun )
  expect_s3_class(out, "tbl_df")
  M = lm( Y0 ~ X1 + X2, data = out )
  expect_equal( as.numeric( coef(M) ),
                c( -100, 100, 0 ), tolerance = 0.1 )
  expect_equal( out$Y0, out$Y1 - 1 )

  # Treatment function
  tx_effect_fun = function(X1, X2) {
    1 + 0.5 * X1 - 0.5 * X2
  }
  out <- gen_df_adv(nc = 20, nt = 10, tx_effect_fun = tx_effect_fun )
  expect_s3_class(out, "tbl_df")
  M = lm( Y1 - Y0 ~ X1 + X2, data = out )
  expect_equal( as.numeric( coef(M) ),
                c( 1, 0.5, -0.5 ), tolerance = 0.1 )


  # heteroskedastic function
  f0_fun <- function( x1, x2 ) {
    0
  }
  f0_sd_fun <- function( x1, x2 ) {
    pmax( 0, 10 + 5 * x1 )
  }
  out <- gen_df_adv(nc = 200, nt = 1000,
                    f0_fun = f0_fun,
                    f0_sd_fun = f0_sd_fun )
  expect_s3_class(out, "tbl_df")
  M = lm( Y0 ~ X1 + X2, data = out )
  coef(M)
  expect_equal( broom::tidy(M)$p.value < 0.001, c( FALSE, FALSE, FALSE ) )
  #expect_equal( as.numeric( coef(M) ),
  #              c( 0, 0, 0 ), tolerance = 0.3 )
  M = lm( noise^2 ~ X1 + X2, data = out )
  expect_equal( broom::tidy(M)$p.value < 0.001, c( TRUE, TRUE, FALSE ) )

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
