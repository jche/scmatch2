

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





test_that( "code at bottom of wrappers just runs", {

  df <- CSM:::gen_one_toy()

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

  expect_true( is.data.frame(res) )

})
