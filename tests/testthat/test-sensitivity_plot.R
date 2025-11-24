

test_that("sensitivity plot does not crash", {

  # Hand-built toy for exploration

  # Make a 2-d covariate example where we scale
  test_df <- tribble( ~X1, ~X2, ~Z, ~Y,
                      1,   3,  1,  1,
                      1,   4,  0,  1,
                      2,   3,  0,  1,

                      11,  0,  1,  0.5,
                      15, 14,  0,  0,

                      4,   2,  1,  1,
                      5,   3,  0,  2,
                      6,   3,  0,  -1,
                      5,   4,  0,  3,

                      6,   12,  1,  1,
                      4,   12,  0,  0,
                      6,   9,  0,  2,

                      3,  13,  1,  4,

                      10, 10, 0, 1,
                      10, 10, 1, 6,

                      10, 15, 0, 1,
                      10, 15, 1, 6

  ) %>%
    group_by( Z ) %>%
    mutate( ID = 1:n() + 100 * Z,
            Y = rnorm( n() ) ) %>%
    ungroup() %>%
    relocate( ID )

  test_df$Wex = 1:nrow(test_df)

  if ( FALSE ) {
    test_df

    ggplot( test_df, aes( X1, X2, col=as.factor(Z), pch=as.factor(Z) ) ) +
      geom_point( size = 3 ) +
      coord_fixed()


  }

  res <- get_cal_matches(df=test_df,
                         Z ~ X1 + X2,
                         rad_method = "fixed",
                         scaling = 1,
                         caliper = 1, warn=FALSE )
  res

  rs = result_table(res, return = "sc_units")
  expect_true( !( "Wex" %in% colnames(rs) ) )

  rs = result_table(res, return="agg_co_units" )
  expect_true(  "Wex" %in% colnames(rs) )
  rs

  if ( FALSE ) {
    result_table(res)


    result_table( res, return = "sc_units" )
    result_table( res, return = "agg_co_units" )


    get_ATT_estimate( res )

  }


  R = 5
  plt <- sensitivity_plot( res, test_df, R = R )

  expect_true( is_ggplot(plt) )

  atts <- attr( plt, "table" )
  expect_true( is.data.frame(atts) )
  expect_equal( nrow(atts), R )

  if ( FALSE ) {

    plt

    # OLD STUFF--Ignore

    ggplot( atts, aes( caliper, ESS_C ) ) +
      geom_line() +
      geom_hline( yintercept = 0, lty=2 )

    aa = result_table( sen[[1]], return = "agg" ) %>%
      filter( Z == 0 )
    aa
    sum( aa$weights )
    ( sum( aa$weights )^2 ) / sum( aa$weights^2 )

    # TODO: 1:1 matching should have ESS = 5, but it is not.  Why?
    # THere is one duplicate unit---shouldn't that get us above 4 though? Or weights drive it lower I guess?
    # But then why the dip at caliper = 2.5?
    #
    #

    result_table( sen[[1]] )
    result_table( sen[[1]], return = "sc_units" )
    result_table( sen[[1]], return = "agg" )


  }


})
