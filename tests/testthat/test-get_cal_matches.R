


test_that("get_cal_matches works well", {

  test_df <-
    data.frame(Z=c(1,0,0,0,1),
               X=c(0,0.5,0.8,3,1.6) )
  test_df
  res <- CSM:::get_cal_matches(data = test_df,
                               covs = "X",
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               scaling = 1)
  res$matches
  expect_equal( nrow( result_table(res) ), 5 )



  test_df <-
    data.frame(Z=c(1,0,0,0,1),
               X1=c(0,0.5,0.8,3,1.6),
               X2=c(0,0,0,0,0),
               GG = rnorm(5) )
  test_df2 = test_df
  test_df2$X2 = 1
  test_df2$Z = 0
  test_df = bind_rows( test_df, test_df2 )
  test_df

  res <- CSM:::get_cal_matches(data = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               scaling = c( 1, 5 ) )
  res$matches
  result_table( res, "agg" )
  expect_equal( nrow( result_table(res, "sc", outcome="GG") ), 4 )

  # And without outcomes
  expect_equal( nrow( result_table(res, "agg") ), 4 )
  expect_equal( nrow( result_table(res, "sc") ), 4 )


  res <- CSM:::get_cal_matches(data = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               scaling = c( 1, 5 ) )
  res$matches
  expect_equal( nrow( result_table( res ) ), 5 )

  # Test the various helper functions
  df = as.data.frame( res )
  df
  expect_equal( nrow( res ), 5  )
  #expect_equal( res[ 1, 4 ], "U1" )
  #expect_equal( res[[ 1, 1]], 1 )
  expect_true( is.data.frame( res[ 1, ] ) )

  expect_output( p <- summary( res, outcome="GG" ), "Subclass sizes")
  # summary(res)

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


  res <- CSM:::get_cal_matches(data = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               est_method = "scm",
                               scaling = 1)
  res$treatment_table
  # Expect subclass all starts with U in treatment table:
  expect_true( all( grepl( "^U", res$treatment_table$subclass ) ) )


  rs <- result_table( res, "agg" )
  rs

  # Expect no subclass in results since units aggregated
  expect_true( all( is.na( rs$subclass ) ) )

  expect_true( all( res$treatment_table$id %in% rs$id ) )


  # different aggregation
  result = result_table( res, "all" )
  expect_true( all( grepl( "^U", res$result$subclass ) ) )
  expect_true( all( !is.na( res$result$dist ) ) )
  expect_equal( nrow( result ), 9 )


  result = result_table( res, "sc" )
  expect_true( all( grepl( "^U", result$subclass ) ) )
  expect_equal( nrow( result ), 8 )



  # Testing switching aggregation ----

  r_co <- agg_co_units( res, outcome="Y" )
  expect_equal( r_co, rs )

  r_syn <- agg_sc_units( res, outcome="Y" )
  expect_equal( r_syn, result_table(res, "sc", outcome="Y") )
  r_avg <- agg_avg_units( res, outcome="Y" )
  r_avg
  expect_equal( r_avg[2,], r_syn[2,] )

})





test_that("adaptive matching gets right units", {

  test_df <-tibble::tribble(
    ~Z,    ~X,
    1,     0,
    0,   0.1,
    0,   0.8,
    0,     3,
    1,     1,
    0,   0.6,
    0,     2,
    1,    10,
    0,  10.1,
    0, 10.05,
    1,    15,
    0, 14.8,
    0, 16,
    0, 16,
    1, 20,
    0, 19.9,
    0, 20.05,
    0, 20.1,
    0, 20.2
  )

  test_df

  scmatches <- get_cal_matches(data=test_df,
                               covs = "X",
                               treatment = "Z",
                               rad_method="adaptive",
                               scaling = 1,
                               k = 2,
                               caliper = 0.2)
  scmatches

  expect_equal( as.numeric( scmatches$adacalipers ),
                c( 0.6, 0.4, 0.2, 1.0, 0.2 ) )

  expect_equal( names(scmatches$adacalipers ),
                scmatches$treatment_table$id )


  expect_equal( scmatches$treatment_table$nc,
                c(2,2,2,3,4) )

  result_table( scmatches )
  result_table( scmatches, return = "sc_units" )

  dt <- get_distance_table( scmatches )
  dt
  expect_true( all( dt$SCM <= dt$closest ) )
  expect_true( all( dt$closest <= dt$average ) )

  expect_equal( dt$SCM,
                c( 0.1, 0.2, 0.05, 0, 0 ) )
})
