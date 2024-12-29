

test_that("get_att_ests works well",{
  test_matched_weighted_df <-
    data.frame(Z= c(1,0,0,0,1),
               Y = c(1,-1,1,0,1),
               weights = c(1,0.5,0.5,1,1) )
  test_matched_weighted_df
  atts <- CSM:::get_att_ests(test_matched_weighted_df)
  atts
  expect_equal( atts, 1 )
})


test_that("get_att_bal should work for a example dataset with either
          Z in 0,1 or Z in F,T, but not in other numerical Z",{

  source( here::here( "scripts/datagen/gen_six_points.R" ) )
  df_six_points <- gen_six_points()
  covs <- c("X1", "X2")
  zform1 <-
    as.formula(paste0("Z ~ ",
                      paste0(covs, collapse="+")))

  test_att_bal = CSM:::get_att_bal(d = df_six_points,
                     form = zform1,
                     tols = rep(0.01, length(covs)))
  expect_equal(test_att_bal, 9.0399,tolerance=1e-3)

  df_Z_logical <-
    df_six_points %>%
    mutate(Z=as.logical(Z))
  test_att_bal = CSM:::get_att_bal(d = df_Z_logical,
                     form = zform1,
                     tols = rep(0.01, length(covs)))
  expect_equal(test_att_bal, 9.0399,tolerance=1e-3)

  # Make sure treatment is binary is checked
  df_Z_random <- df_six_points %>%
    ungroup() %>%
    mutate( Z = as.numeric(1:6))
  expect_error(
    CSM:::get_att_bal(d = df_Z_random,
                form = zform1,
                tols = rep(0.01, length(covs))),
    "Treatment variable must be either 0, 1 or FALSE, TRUE")
})



test_that( "cem works", {
  df <- CSM:::gen_one_toy()

  cem <- CSM:::get_att_cem(df, num_bins=5, est_method="scm")

  expect_true( is.numeric(cem) )
})

test_that("get_matches works with different matching types", {
  test_df <- data.frame(
    Z = c(1, 0, 0, 0, 1),
    X = c(0, 0.5, 0.8, 3, 1.6)
  )
  scaling <- 1  # Example scaling factor

  res_toy <- get_matches(
    matching_type = "maximum_fixed_scm",
    df_dgp = test_df,
    scaling = scaling
  )
  expect_equal(nrow(res_toy), 5)  # Ensure 4 rows in the result

  # Test data for otsu
  test_df_otsu <- data.frame(
    Z = c(1, 0, 0, 0, 1),
    X1 = c(0, 0.5, 0.8, 3, 1.6),
    X2 = c(0, 0, 0, 0, 0)
  )
  test_df_otsu <- bind_rows(
    test_df_otsu,
    transform(test_df_otsu, X2 = 1, Z = 0)
  )

  # Test for 'euclidean_knn'
  res_otsu <- get_matches(
    matching_type = "euclidean_knn",
    df_dgp = test_df_otsu,
    scaling = scaling
  )
  expect_equal(nrow(res_otsu), 18)

  # Test for invalid matching_type
  expect_error(
    get_matches(
      matching_type = "invalid_type",
      df_dgp = test_df,
      scaling = scaling
    ),
    "Invalid matching_type"
  )
})



test_that( "all other methods run and give ATT estimates", {

  set.seed( 40440 )
  df <- CSM:::gen_one_toy()
  df

  df %>%
    ggplot(aes(X1,X2)) +
    geom_point(aes(pch=as.factor(Z), color=Y)) +
    scale_color_continuous(low="orange", high="blue") +
    theme_classic() +
    labs(pch = "Treated",
         color = latex2exp::TeX("$f_0(X)$"),
         x = latex2exp::TeX("$X_1$"),
         y = latex2exp::TeX("$X_2$")) +
    facet_wrap(~Z)

  df %>%
    filter(Z) %>%
    summarize(eff = mean(Y1-Y0))

  res <- CSM:::run_all_methods( df )
  res
  expect_true( is.data.frame(res) )


  df = filter( df, Z==1 )
  df2 = df %>%
    mutate( Z = 0,
            Y = Y - 1 )

  df = bind_rows(df, df2)
  df

  # Looking at CSM local averaging issue
  preds_csm <- get_cal_matches(
    df = df,
    metric = "maximum",
    scaling = 100,
    #caliper = 0.001,
    rad_method = "adaptive",
    est_method = "scm",
    return = "sc_units")
  preds_csm$treatment_table
  head( preds_csm$matches )

  tt <- CSM:::get_att_csm(df, scaling=100)
  expect_equal( tt, 1 )

  # NOTE: CSM is biased here due to radial matching keeping lots of
  # units due to large radius
  res <- CSM:::run_all_methods( df, extrapolate = FALSE )
  res <- filter( res, !is.na(att) )
  res <- filter( res, method != "csm" )
  res <- filter( res, method != "cem" ) # It extrapolates and messes up?
  expect_equal( res$att, rep( 1, nrow(res) ), tolerance=1e-2 )
})
