
## Main test: gen_matches


test_that("gen_matches returns correct matches
          under fixed caliper and two treatments with
          one shared control", {

            test_df <-
              data.frame(Z=c(1,0,0,0,1),
                         X=c(0,0.5,0.8,3,1.6))
  test_df
            res <- gen_matches(data=test_df,
                               covs = "X",
                               treatment = "Z",
                               caliper = 1)

            res$matches

            res_matched_pair_1 <- as.matrix(res$matches[[1]])
            expected_matched_pair_1 <-
              data.frame(
                id = c( "tx1", "co1", "co2" ),
                Z = c(1, 0, 0),
                X = c(0, 0.5, 0.8),
                dist = c(0, 0.5, 0.8),
                subclass = c("tx1", "tx1", "tx1")
              )
            #rownames(expected_matched_pair_1) <- c('1','2','3')
            expect_equal(res_matched_pair_1,
                         as.matrix(expected_matched_pair_1))

            res_matched_pair_2 <- as.matrix(res$matches[[2]])
            res_matched_pair_2
            expected_matched_pair_2 <-
              data.frame(
                id = c("tx2", "co2" ),
                Z = c(1, 0),
                X = c(1.6, 0.8),
                dist = c(0, 0.8),
                subclass = c("tx2", "tx2")
              )
            #rownames(expected_matched_pair_2) <- c('1','3')
            expect_equal(res_matched_pair_2,
                         as.matrix(expected_matched_pair_2))
          })



test_that("gen_matches did well in different scaling",{

  test_df <-
    data.frame(Z=c(1,0,0),
               X=c(0,0.5,0.8))
  res <- gen_matches(data=test_df,
                     covs = "X",
                     treatment = "Z",
                     scaling = 1.5,
                     caliper = 1)

  res_matched_pair_1 <- as.matrix(res$matches[[1]])
  res_matched_pair_1
  expected_matched_pair_1 <-
    data.frame(
      id = c("tx1", "co1"),
      Z = c(1, 0),
      X = c(0, 0.5), # X retains in the original scale
      dist = c(0, 0.75), # but distance is computed in scaled X
      subclass = c("tx1", "tx1")
    )
 # rownames(expected_matched_pair_1) <- c('1','2')
  expect_equal(res_matched_pair_1,
               as.matrix(expected_matched_pair_1))

  # The above should be equivalent to scaling = 1, but caliper = 1/1.5
  res <- gen_matches(data=test_df,
                     covs = "X",
                     treatment = "Z",
                     scaling = 1,
                     caliper = 1/1.5)

  res_matched_pair_1 <- as.matrix(res$matches[[1]])
  expected_matched_pair_1 <-
    data.frame(
      id = c("tx1", "co1"),
      Z = c(1, 0),
      X = c(0, 0.5),
      dist = c(0, 0.5), # distance is retained in original scale, but the last unit got unmatched
      subclass = c("tx1", "tx1")
    )
  #rownames(expected_matched_pair_1) <- c('1','2')
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


  res <- gen_matches(data=test_df,
                     covs = c( "X1", "X2" ),
                     treatment = "Z",
                     scaling = 1,
                     caliper = 1)
  res
  expect_equal( as.numeric(res$adacalipers), c( 1,5,1,2 ) )


  res <- gen_matches(data=test_df,
                     covs = c( "X1", "X2" ),
                     treatment = "Z",
                     metric = "euclidean",
                     scaling = 1,
                     caliper = 1)
  res$matches
  expect_equal( as.numeric(res$adacalipers), c( 1,sqrt(5^2+3^2),sqrt(2) ,2 ) )

  res <- get_cal_matches(data=test_df,
                         covs = c( "X1", "X2" ),
                         treatment = "Z",
                         metric = "euclidean",
                         scaling = 1,
                         caliper = 5)
  rt = result_table( res, "sc", outcome = "Y" )
  expect_equal( nrow( rt ), 8 )


  resf <- get_cal_matches(data=test_df,
                          covs = c( "X1", "X2" ),
                          treatment = "Z",
                          metric = "euclidean",
                          rad_method = "fixed",
                          scaling = 1,
                          caliper = 5, warn=FALSE)

  rtf = result_table(resf, "sc", outcome = "Y" )
  expect_equal( nrow( rtf ), 6 )

  rtf
  rt
  expect_equal( rtf[5:6,]$Y,
                rt[7:8,]$Y, tolerance=1e-5 )



  res <- get_cal_matches(data=test_df,
                         covs = c( "X1", "X2" ),
                         treatment = "Z",
                         scaling = c(1, 0.25),
                         rad_method="fixed",
                         caliper = 1, warn = FALSE)
  res
  expect_equal( as.numeric(res$adacalipers), c( 1, NA, 1, 1 ) )


  # Scaling and caliper gives same thing?
  retx2 <- get_cal_matches(data=test_df,
                           covs = c( "X1", "X2" ),
                           treatment = "Z",
                           scaling = c(4, 1),
                           rad_method = "fixed",
                           caliper = 1/0.25, warn=FALSE)


  expect_equal( result_table( res )$weights, result_table( retx2 )$weights )
  expect_equal( result_table( res )$ID, result_table( retx2 )$ID )
})





test_that("gen_match adaptive caliper", {

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


  res <- gen_matches(data = test_df,
                               covs = c( "X1", "X2" ),
                               treatment = "Z",
                               metric = "maximum",
                               caliper = 1,
                               rad_method = "adaptive",
                               scaling = 1)

  treatment_table = CSM:::make_treatment_table(res)
  expect_equal( as.numeric(treatment_table$adacal),
                c( 1, 5, 1, 2 ) )
  expect_equal( as.numeric(treatment_table$max_dist),
                c( 0.5, 5, 1, 2 ) )
  res


})




test_that("get_radius_size should get the correct output",{
  caliper=1
  dm = matrix(c(1.5, 2, 1.6, 1.5,
                  1.3, 1.9, 1.8, 0.7,
                  0.9, 0.5, 0.3, 0.2),
                byrow = TRUE, nrow=3)
  dm
  res <- get_radius_size(dm=dm,
                         rad_method="adaptive",
                         caliper=caliper)
  res
  expected_res <- c(1.5, 1, 1)
  expect_equal(res, expected_res)

  res <- get_radius_size(dm=dm,
                         rad_method="knn")
  res
  expect_equal( res, apply(dm,1,min) )

  res2 <- get_radius_size(dm=dm,
                          rad_method="1nn")
  expect_equal( res2, res )


  res2 <- get_radius_size(dm=dm,
                          rad_method="knn-capped", caliper = 1)
  res2
  expect_equal( res2, pmin( apply(dm,1,min), 1 ) )


  res <- get_radius_size(dm=dm,
                         rad_method="knn", k = 2)
  res
  expect_equal( res, apply(dm,1,function(x) sort(x)[2] ) )

  res <- get_radius_size(dm=dm,
                         rad_method="targeted",
                         caliper = 1, k = 2)
  res
  expect_equal( res, c( 1.5, 1.0, 0.3 ) )



})




test_that("set_NA_to_unmatched_co should get the correct output",{
  test_df <-
    data.frame(Z=c(1,1, 0,0,0),
               X=c(0,3.2, 0.5,0.8,3))
  dm_uncapped <- CSM:::gen_dm(data=test_df,
                              covs="X",
                              treatment="Z",
                              scaling=1,
                              metric="maximum")
  dm_uncapped
  radius_sizes <- c(1, 0.1)
  res <- CSM:::set_NA_to_unmatched_co(dm_uncapped, radius_sizes)
  res
  expected_res <-  t(matrix(c(0.5, 0.8, NA),ncol=1))
  colnames(expected_res) <- c(3,4,5)
  rownames(expected_res) <- 1

  attr( res, "covariates" ) <- NULL
  expect_equal(res[1,, drop=FALSE], expected_res)
  expect_equal( as.numeric(res[2,]),
                as.double(c(NA,NA,NA) ) )

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
  res <- CSM:::get_matched_co_from_dm_trimmed(data=test_df,
                                              dm_trimmed=dm_trimmed,
                                              treatment=treatment)
  # Convert the actual result (list element) to a standard data frame
  actual_df <- as.data.frame(res[[1]])

  # Ensure row names are removed/reset on both for a clean comparison
  rownames(actual_df) <- NULL
  rownames(expected_res) <- NULL

  expect_equal(actual_df, expected_res)
})

test_that("get_matched_co_from_dm_trimmed handles duplicate covariates (bootstrap ties) correctly", {

  # 1. Setup: Create a dataset where control 2 and control 3 are identical in X
  # but have unique IDs (simulating a bootstrap draw with re-id)
  test_df <- tibble(
    id = c("tx1", "co1", "co2"),
    X  = c(0.5, 0.6, 0.6), # Tied covariate values
    W  = c(1, 0, 0),       # Treatment indicator
    Y  = c(10, 11, 12)     # Distinct outcomes
  )

  # 2. Setup: Distance matrix where Treated 1 matches both controls
  # dm_trimmed rows = treated, columns = controls
  dm_trimmed <- matrix(c(0.1, 0.1), nrow = 1)
  # In your CSM logic, column names often refer to the original control indices
  colnames(dm_trimmed) <- c("2", "3")

  # 3. Execution
  res <- get_matched_co_from_dm_trimmed(test_df, dm_trimmed, "W")

  # 4. Assertions
  # Check that the list returned is valid
  expect_type(res, "list")
  expect_length(res, 1)

  match_table <- res[[1]]

  # A. Check for NAs: The core issue you reported
  expect_false(any(is.na(match_table$id)),
               info = "IDs should not be NA. If NA, row name indexing is likely broken.")
  expect_false(any(is.na(match_table$X)),
               info = "Covariates should not be NA.")

  # B. Check structure: 1 Treated + 2 Controls = 3 Rows
  expect_equal(nrow(match_table), 3)

  # C. Check Content: Ensure Control 'co2' (row 3 of match_table) has the correct Y
  # If indexing is broken, row 3 might accidentally be a copy of row 2
  expect_equal(match_table$id, c("tx1", "co1", "co2"))
  expect_equal(match_table$Y, c(10, 11, 12))

  # D. Check subclass assignment
  expect_true(all(match_table$subclass == "tx1"))
})


test_that("est_weights works fine", {
  test_df <-
    data.frame(Z=c(1,0,0,0, 1,0,0),
               X=c(0,0.5,0.8,3, 1,0.6,2))
  test_df

  scmatches <- gen_matches(data=test_df,
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

  scweights2 <- est_weights(scmatches, est_method = "average")$matches
  expect_equal( length( scweights2 ), 2 )
  expect_equal( scweights2[[2]]$weights, c( 1, rep(0.25, 4 ) ) )

  scweights2 <- est_weights(scmatches, est_method = "scm")
  expect_equal( scweights, scweights2$matches )

})






