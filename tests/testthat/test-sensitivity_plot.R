
# Hand-built toy for exploration

# Make a 2-d covariate example where we scale
test_df <- tribble( ~WX1, ~WX2, ~ZZ, ~YY,
                    1,   3,  1,  1,
                    1,   4,  0,  2,
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
  group_by( ZZ ) %>%
  mutate( ID = 1:n() + 100 * ZZ,
          Y2 = rnorm( n() ) ) %>%
  ungroup() %>%
  relocate( ID )

test_df$Wex = 1:nrow(test_df)


if ( FALSE ) {
  test_df

  ggplot( test_df, aes( X1, X2, col=as.factor(Z), pch=as.factor(Z) ) ) +
    geom_point( size = 3 ) +
    coord_fixed()


}




test_that("sensitivity plot does not crash", {


  res <- get_cal_matches(data=test_df,
                         ZZ ~ WX1 + WX2,
                         rad_method = "fixed",
                         scaling = 1,
                         caliper = 1, warn=FALSE )
  res

  rs = result_table(res, return = "sc_units", outcome="YY")
  expect_true( !( "Wex" %in% colnames(rs) ) )

  rs = result_table(res, return="agg_co_units" )
  expect_true(  "Wex" %in% colnames(rs) )
  rs

  if ( FALSE ) {
    result_table(res)
    result_table( res, return = "sc_units" )
    result_table( res, return = "agg_co_units" )
    estimate_ATT( res )
  }

  res2 <- update_matches( res, test_df, rad_method = "adaptive" )
  res2
  pp <- sensitivity_table( res2, outcome="YY" )
  pp
  expect_true( is.data.frame(pp) )
  expect_true( all( c( "ATT", "SE", "ESS_C" ) %in% colnames(pp) ) )


  aa <- update_matches( res, test_df, est_method = "average", warn=FALSE )
  ee <- estimate_ATT( aa, outcome="YY" )
  pa = filter(pp, Estimate == "FATT_raw" )
  pa
  ee
  expect_equal( dplyr::select( ee, ATT, SE, N_T, ESS_C ),
                dplyr::select( pa, ATT, SE, N_T, ESS_C ) )



  aa <- update_matches( res, test_df, rad_method = "1nn" )
  ee <- estimate_ATT( aa, outcome="YY" )
  pa = filter(pp, Estimate == "ATT_1nn" )
  pa
  ee
  expect_equal( dplyr::select( ee, ATT, SE, N_T, ESS_C ),
                dplyr::select( pa, ATT, SE, N_T, ESS_C ) )


} )


test_that("caliper sensitivity plot works", {

  R = 5
  res <- get_cal_matches(data=test_df,
                         ZZ ~ WX1,
                         rad_method = "adaptive",
                         scaling = 2,
                         caliper = 0.2, warn=FALSE )

  tbl <- caliper_sensitivity_table( res, outcome="YY", test_df, R = R, min_cal = 0, max_cal = 3 )
  expect_true( is.data.frame( tbl ) )


  plt <- caliper_sensitivity_plot( tbl, outcome="YY", test_df, R = R, min_cal = 0, max_cal = 3 )
  if ( FALSE ) {
    plt
  }
  expect_true( is_ggplot(plt) )

  plt2 <- caliper_sensitivity_plot( res, outcome="YY", test_df, R = R, min_cal = 0, max_cal = 3 )
  expect_equal( attr( plt, "table" ),
                attr( plt2, "table" ) )


  plt3 <- caliper_sensitivity_plot_stats(plt2)
  plt3

  plt3 <- caliper_sensitivity_plot_stats(plt2,
                                         vars = c("ATT", "SE", "ESS_C" ) )
  plt3

  atts <- attr( plt3, "table" )
  atts
  expect_true( is.data.frame(atts) )
  expect_equal( nrow(atts), R * 6 )

  # SHould have full ribbon on plot
  plt2 <- caliper_sensitivity_plot( atts, outcome="YY", focus="ATT", test_df, R = R, min_cal = 0, max_cal = 3 )
  plt2
  expect_true( is_ggplot(plt2) )


  # can facet with multiple foci
  plt2 <- caliper_sensitivity_plot( atts, outcome="YY", focus=c("ATT", "FATT"), test_df, R = R, min_cal = 0, max_cal = 3 ) +
    facet_wrap( ~ feasible )
  plt2
  expect_true( is_ggplot(plt2) )

})
