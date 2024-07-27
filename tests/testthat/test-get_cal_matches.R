


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
  res$matches
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
  res$matches
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
  res$matches
  expect_equal( nrow( res$result ), 5 )


})




test_that("get_cal_matches adaptive caliper", {

  test_df <- tribble( ~X1, ~X2, ~Z, ~Y,
                      1,   3,  1,  1,
                      1,   3.5, 0, 1,
                      1.5, 3.5, 0, 1,

                      11,  0,  1,  0,
                      15, 14,  0,  0,

                      4,   2,  1,  1,
                      5,   3,  0,  0,
                      6,   3,  0,  0,
                      5,   4,  0,  0,

                      6,   12,  1,  1,
                      4,   12,  0,  0,
                      6,   9,  0,  2) %>%
    mutate( ID = 1:n() + 100 * Z,
            Y = rnorm( n() ) ) %>%
    relocate( ID )


  res <- CSM:::get_cal_matches(df = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               return = "agg_co_units",
                               scaling = 1)
  res$treatment_table
  expect_equal( res$treatment_table$adacal, c( 1, 5, 1, 2 ) )
  expect_equal( res$treatment_table$max_dist, c( 0.5, 5, 1, 2 ) )
  res


})







test_that("id and subclass ID make sense", {

  test_df <- tribble( ~X1, ~X2, ~Z, ~Y,
                      1,   3,  1,  1,
                      1,   3.5, 0, 1,
                      1.5, 3.5, 0, 1,

                      11,  0,  1,  0,
                      15, 14,  0,  0,

                      4,   2,  1,  1,
                      5,   3,  0,  0,
                      6,   3,  0,  0,
                      5,   4,  0,  0,

                      6,   12,  1,  1,
                      4,   12,  0,  0,
                      6,   9,  0,  2) %>%
    mutate( ID = 1:n() + 100 * Z,
            Y = rnorm( n() ) ) %>%
    relocate( ID )


  res <- CSM:::get_cal_matches(df = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               return = "agg_co_units",
                               scaling = 1)
  res$treatment_table
  # Expect subclass all starts with S in treatment table:
  expect_true( all( grepl( "^U", res$treatment_table$subclass ) ) )


  rs <- result_table( res )
  rs
  # Expect no subclass in results since units aggregated
  expect_true( all( is.na( rs$subclass ) ) )

  expect_true( all( res$treatment_table$id %in% rs$id ) )


  # different aggregation
  res <- CSM:::get_cal_matches(df = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               return = "all",
                               #return = c("sc_units", "agg_co_units", "all"),

                               scaling = 1)
  res$result
  expect_true( all( grepl( "^U", res$result$subclass ) ) )
  expect_true( all( !is.na( res$result$dist ) ) )
  expect_equal( nrow( res$result ), 9 )


  res_sc_units <- CSM:::get_cal_matches(df = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               return = "sc_units",
                               #return = c("sc_units", "agg_co_units", "all"),

                               scaling = 1)
  res_sc_units$result
  expect_true( all( grepl( "^U", res_sc_units$result$subclass ) ) )
  expect_equal( nrow( res_sc_units$result ), 8 )



  # Testing switching aggregation ----

  r_co <- agg_co_units( res )
  expect_equal( r_co, rs )
  r_syn <- agg_sc_units( res )
  expect_equal( r_syn, result_table(res_sc_units) )
  r_avg <- agg_avg_units( res )
  r_avg
  expect_equal( r_avg[2,], r_syn[2,] )
  expect_equal( res_sc_units$result, r_syn )
})


