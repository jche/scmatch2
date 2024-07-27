


test_that("get_cal_matches works well", {
  test_df <-
    data.frame(Z=c(1,0,0,0,1),
               X=c(0,0.5,0.8,3,1.6) )
  res <- CSM:::get_cal_matches(df = test_df,
                               covs = "X",
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               return = "agg_co_units",
                               scaling = 1)
  res
  expect_equal( nrow( res$result ), 4 )



  test_df <-
    data.frame(Z=c(1,0,0,0,1),
               X1=c(0,0.5,0.8,3,1.6),
               X2=c(0,0,0,0,0) )
  test_df2 = test_df
  test_df2$X2 = 1
  test_df2$Z = 0
  test_df = bind_rows( test_df, test_df2 )
  test_df

  res <- CSM:::get_cal_matches(df = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               return = "agg_co_units",
                               scaling = c( 1, 5 ) )
  res
  expect_equal( nrow( res$result ), 4 )


  res <- CSM:::get_cal_matches(df = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               return = "all",
                               scaling = c( 1, 5 ) )
  res
  expect_equal( nrow( res$result ), 5 )

})



