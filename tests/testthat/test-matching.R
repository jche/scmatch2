
## Main test: gen_matches


test_that("gen_matches returns correct matches
          under fixed caliper and two treatments with
          one shared control", {

            test_df <-
              data.frame(Z=c(1,0,0,0,1),
                         X=c(0,0.5,0.8,3,1.6))

            res <- gen_matches(df=test_df,
                               covs = "X",
                               treatment = "Z",
                               caliper = 1)

            res$matches
            res_matched_pair_1 <- as.matrix(res$matches[[1]])
            expected_matched_pair_1 <-
              data.frame(
                Z = c(1, 0, 0),
                X = c(0, 0.5, 0.8),
                dist = c(0, 0.5, 0.8),
                subclass = c(1, 1, 1)
              )
            rownames(expected_matched_pair_1) <- c('1','2','3')
            expect_equal(res_matched_pair_1,
                         as.matrix(expected_matched_pair_1))

            res_matched_pair_2 <- as.matrix(res$matches[[2]])
            expected_matched_pair_2 <-
              data.frame(
                Z = c(1, 0),
                X = c(1.6, 0.8),
                dist = c(0, 0.8),
                subclass = c(2, 2)
              )
            rownames(expected_matched_pair_2) <- c('1','3')
            expect_equal(res_matched_pair_2,
                         as.matrix(expected_matched_pair_2))
          })



test_that("gen_matches did well in different scaling",{
  test_df <-
    data.frame(Z=c(1,0,0),
               X=c(0,0.5,0.8))
  res <- gen_matches(df=test_df,
                     covs = "X",
                     treatment = "Z",
                     scaling = 1.5,
                     caliper = 1)

  res_matched_pair_1 <- as.matrix(res$matches[[1]])
  expected_matched_pair_1 <-
    data.frame(
      Z = c(1, 0),
      X = c(0, 0.5), # X retains in the original scale
      dist = c(0, 0.75), # but distance is computed in scaled X
      subclass = c(1, 1)
    )
  rownames(expected_matched_pair_1) <- c('1','2')
  expect_equal(res_matched_pair_1,
               as.matrix(expected_matched_pair_1))

  # The above should be equivalent to scaling = 1, but caliper = 1/1.5
  res <- gen_matches(df=test_df,
                     covs = "X",
                     treatment = "Z",
                     scaling = 1,
                     caliper = 1/1.5)

  res_matched_pair_1 <- as.matrix(res$matches[[1]])
  expected_matched_pair_1 <-
    data.frame(
      Z = c(1, 0),
      X = c(0, 0.5),
      dist = c(0, 0.5), # distance is retained in original scale, but the last unit got unmatched
      subclass = c(1, 1)
    )
  rownames(expected_matched_pair_1) <- c('1','2')
  expect_equal(res_matched_pair_1,
               as.matrix(expected_matched_pair_1))
})



test_that("gen_matches do the correct thing with non-uniform scaling", {
  # Make a 2-d covariate example where we scale
  test_df <- tribble( ~X1, ~X2, ~Z, ~Y,
                      1,   3,  1,  1,
                      1,   4,  0,  1,
                      2,   3,  0,  1,

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


  res <- gen_matches(df=test_df,
                     covs = c( "X1", "X2" ),
                     treatment = "Z",
                     scaling = 1,
                     caliper = 1)
  res
  expect_equal( res$adacalipers, c( 1,5,1,2 ) )


  res <- gen_matches(df=test_df,
                     covs = c( "X1", "X2" ),
                     treatment = "Z",
                     metric = "euclidean",
                     scaling = 1,
                     caliper = 1)
  res$matches
  expect_equal( res$adacalipers, c( 1,sqrt(5^2+3^2),sqrt(2) ,2 ) )

  res <- get_cal_matches(df=test_df,
                         covs = c( "X1", "X2" ),
                         treatment = "Z",
                         metric = "euclidean",
                         scaling = 1,
                         caliper = 5)
  expect_equal( nrow( res$result ), 8 )


  resf <- get_cal_matches(df=test_df,
                          covs = c( "X1", "X2" ),
                          treatment = "Z",
                          metric = "euclidean",
                          rad_method = "fixed",
                          scaling = 1,
                          caliper = 5, warn=FALSE)
  expect_equal( nrow( resf$result ), 6 )

  expect_equal( resf$result[5:6,]$Y,
                res$result[7:8,]$Y )



  res <- get_cal_matches(df=test_df,
                         covs = c( "X1", "X2" ),
                         treatment = "Z",
                         scaling = c(1, 0.25),
                         rad_method="fixed",
                         return = "all",
                         caliper = 1)
  res
  expect_equal( res$adacalipers, c( 1, NA, 1, 1 ) )


  # Scaling and caliper gives same thing?
  res2 <- get_cal_matches(df=test_df,
                          covs = c( "X1", "X2" ),
                          treatment = "Z",
                          scaling = c(4, 1),
                          rad_method = "fixed",
                          return = "all",
                          caliper = 1/0.25)


  expect_equal( res$result$weights, res2$result$weights )
  expect_equal( res$result$ID, res2$result$ID )
})





test_that("get_radius_size should get the correct output",{
  caliper=1
  dm = t(matrix(c(1.5, 2, # should have radius_size 1.5
                  1.3, 0.9, # should have radius_size 1
                  0.9,0.5), # should have radius_size 1
                nrow=2) )
  res <- get_radius_size(dm=dm,
                         rad_method="adaptive",
                         caliper=caliper)
  expected_res <- c(1.5, 1, 1)
  expect_equal(res, expected_res)
})




test_that("set_NA_to_unmatched_co should get the correct output",{
  test_df <-
    data.frame(Z=c(1,0,0,0),
               X=c(0,0.5,0.8,3))
  dm_uncapped <- CSM:::gen_dm(df=test_df,
                              covs="X",
                              treatment="Z",
                              scaling=1,
                              metric="maximum")
  radius_sizes <- c(1)
  res <- CSM:::set_NA_to_unmatched_co(dm_uncapped, radius_sizes)

  expected_res <-  t(matrix(c(0.5, 0.8, NA),ncol=1))
  colnames(expected_res) <- c(2,3,4)
  rownames(expected_res) <- 1

  expect_equal(res, expected_res)
})



test_that("get_matched_co_from_dm should get the correct output",{
  test_df <-
    data.frame(Z=c(1,0,0,0),
               X=c(0,0.5,0.8,3))
  dm_trimmed <-  t(matrix(c(0.5, 0.8, NA),ncol=1))
  colnames(dm_trimmed) <- c(2,3,4)
  rownames(dm_trimmed) <- 1
  treatment = "Z"
  test_df_trt <- test_df %>%
    filter(.data[[treatment]] == 1)
  res <- CSM:::get_matched_co_from_dm_trimmed(df=test_df,
                                              dm_trimmed=dm_trimmed,
                                              treatment=treatment)
  res_values <- as.matrix(res[[1]])
  expected_res <-
    data.frame(
      Z = c(1, 0, 0),
      X = c(0, 0.5, 0.8),
      dist = c(0, 0.5, 0.8),
      subclass = c(1, 1, 1)
    )
  rownames(expected_res) <- c('1','2','3')
  expect_equal(res_values, as.matrix(expected_res))
})




test_that("est_weights works fine", {
  test_df <-
    data.frame(Z=c(1,0,0,0, 1,0,0),
               X=c(0,0.5,0.8,3, 1,0.6,2))
  test_df

  scmatches <- gen_matches(df=test_df,
                           covs = "X",
                           treatment = "Z",
                           caliper = 1)
  scmatches
  expect_equal( length( scmatches$matches ), 2 )

  covs = "X"
  matched_gps = scmatches$matches
  scaling=1
  metric = "maximum"

  scweights <-
    est_weights(covs = covs,
                matched_gps = matched_gps,
                scaling = scaling,
                est_method = "scm",
                metric = metric)

  scweights
  expect_equal( length( scweights ), 2 )
})



