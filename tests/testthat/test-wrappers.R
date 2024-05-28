
test_that("fit_CSM works well",{
  test_df <-
    data.frame(Z=c(1,0,0,0,1),
               X=c(0,0.5,0.8,3,1.6),
               Y = c(1,0,0,0,1))
  res <- fit_CSM(df = test_df,
                         covs = "X",
                         treatment = "Z",
                         metric = "maximum",
                         caliper = 1,
                         rad_method = "adaptive",
                         est_method = "scm",
                         return = "agg_co_units",
                         dist_scaling = 1)
  expect_equal(res, 1)
})


test_that("get_cal_matches works well", {
  test_df <-
    data.frame(Z=c(1,0,0,0,1),
               X=c(0,-0.5,0.8,3,1.6))
  attr(res, "feasible_units")
  res <- get_cal_matches(df = test_df,
                         covs = "X",
                         treatment = "Z",
                         metric = "maximum",
                         caliper = 1,
                         rad_method = "adaptive",
                         est_method = "scm",
                         return = "agg_co_units",
                         dist_scaling = 1)
})

test_that("get_att_ests works well",{
  test_matched_weighted_df <-
    data.frame(Z= c(1,0,0,0,1),
               Y = c(1,0,0,0,1),
               weights = )
  get_att_ests(res)
})


test_that("get_att_bal should work for a example dataset with either
          Z in 0,1 or Z in F,T, but not in other numerical Z",{
  source("./scripts/datagen/gen_six_points.R")
  df_six_points <- gen_six_points()
  covs <- c("X1", "X2")
  zform1 <-
    as.formula(paste0("Z ~ ",
                      paste0(covs, collapse="+")))

  test_att_bal = get_att_bal(d = df_six_points,
                     form = zform1,
                     tols = rep(0.01, length(covs)))
  expect_equal(test_att_bal, 9.0399,tolerance=1e-3)

  df_Z_logical <-
    df_six_points %>%
    mutate(Z=as.logical(Z))
  test_att_bal = get_att_bal(d = df_Z_logical,
                     form = zform1,
                     tols = rep(0.01, length(covs)))
  expect_equal(test_att_bal, 9.0399,tolerance=1e-3)

  df_Z_random <- df_six_points %>%
    ungroup() %>%
    mutate(Z= as.numeric(1:6))
  expect_error(
    get_att_bal(d = df_Z_random,
                form = zform1,
                tols = rep(0.01, length(covs))),
    "Treatment variable must be either 0, 1 or FALSE, TRUE")
})


test_that("get_SL_pred returns correct output structure", {
  mock_SL_fit <- create_mock_SL_fit()
  mock_df_test <- create_mock_df_test()
  X_names <- c("X1","X2")
  result <- get_SL_pred(mock_SL_fit,
                        mock_df_test,
                        X_names)

  expect_true(is.matrix(result) == T)
})

test_that("get_SL_fit returns a correct output data type", {
  result = create_mock_SL_fit()
  expect_true("SuperLearner" %in% class(result))
})
