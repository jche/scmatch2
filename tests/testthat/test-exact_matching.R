

test_that("exact matching defaults work", {

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
                      10, 15, 1, 6,

                      20, 20, 1, 5,
                      20, 20.5, 0, 2

  ) %>%
    group_by( Z ) %>%
    mutate( ID = 1:n() + 100 * Z,
            Y = rnorm( n() ) ) %>%
    ungroup() %>%
    relocate( ID )

  test_df


  res <- get_cal_matches(data=test_df,
                         Z ~ X1 + X2,
                         rad_method = "fixed",
                         scaling = 1,
                         caliper = 0, warn=FALSE )
  rs <- result_table(res)
  rs
  expect_equal( nrow( rs ), 4 )
  expect_equal( rs$dist, c(0,0,0,0) )

  expect_output( summary( res ) )

  res <- get_cal_matches(data=test_df,
                         Z ~ X1 + X2,
                         rad_method = "fixed",
                         scaling = 1,
                         caliper = 1, warn=FALSE )
  rs <- result_table(res)
  rs
  rs2 <- result_table( res, return = "exact" )
  rs2

  expect_equal( nrow( rs2 ), 4 )
  expect_equal( as.numeric(rs2$dist), c(0,0,0,0) )

  expect_output( summary( res ) )


  # And no exact matching?
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

                      3,  13,  1,  4
  ) %>%
    group_by( Z ) %>%
    mutate( ID = 1:n() + 100 * Z,
            Y = rnorm( n() ) ) %>%
    ungroup() %>%
    relocate( ID )

  res <- get_cal_matches(data=test_df,
                         Z ~ X1 + X2,
                         rad_method = "fixed",
                         scaling = 1,
                         caliper = 0, warn=FALSE )
  rs <- result_table(res)
  rs
  expect_equal( nrow( rs ), 0 )

  expect_output( summary( res ) )


  res <- get_cal_matches(data=test_df,
                         Z ~ X1 + X2,
                         scaling = 1,
                         caliper = 0, warn=FALSE )

  tb = result_table( res, "exact" )
  expect_equal( nrow(tb), 0 )




})
