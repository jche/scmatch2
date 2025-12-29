# tests/testthat/test-csm_matches_object.R

# Test the various functionalities of the csm_matches_object


# The various methods for the csm_matches object




test_that("basic object functionality", {

  # Generate example data
  set.seed(4044440)
  dat <- gen_one_toy(nt = 5)

  # Perform matching
  mtch <- get_cal_matches(dat,
                          treatment = "Z",
                          metric = "maximum",
                          scaling = c(1/0.2, 1/0.2),
                          caliper = 0.25,
                          rad_method = "adaptive",
                          est_method = "scm")

  # View matching results
  mtch

  expect_true( is.csm_matches( mtch ) )
  expect_true( !is.csm_matches( as.data.frame( mtch ) ) )


  a = mtch[ 1, ]
  expect_true( is.data.frame( a ) )

  expect_true( mtch[[ 1, 4 ]] )

  dim(mtch)
  expect_true( is.numeric( dim( mtch )[1] ) )
  expect_true( is.numeric( dim( mtch )[2] ) )

  dd = as.data.frame( mtch )
  expect_true( is.data.frame( dd ) )
  expect_equal( nrow( dd ), nrow( result_table( mtch ) ) )

  expect_output( ss <- print( mtch ) )

  expect_output( ss <- summary( mtch ) )


  tb <- caliper_table( mtch )
  expect_true( is.data.frame( tb ) )
  expect_true( all( c("id","subclass","max_dist","adacal") %in% colnames( tb ) ) )
  expect_true( all( tb$id %in% mtch$treatment_table$id ) )

  ff = feasible_units( mtch )
  ff
  expect_true( all( ff$max_dist <= 0.25 ) )

  dl = tb %>%
    dplyr::filter( !( id %in% feasible_units( mtch )$id ) )
  dl
  expect_true( all( dl$max_dist > 0.25 ) )


  fs = feasible_unit_subclass( mtch )
  expect_true( all( fs %in% ff$id ) )

  um <- unmatched_units( mtch )
  expect_true( is.data.frame( um ) )
  expect_equal( nrow( um ), 0 )

  expect_warning(
    umt <- update_matches( mtch, dat, caliper = 0.35,
                           rad_method = "fixed",
                           est_method = "average" )
  )
  umt
  expect_true( is.csm_matches( umt ) )

  rs = result_table( umt )
  rs
  expect_true( sum( rs$Z ) < nrow( mtch$treatment_table ) )


  rs = result_table( umt, return = "exact" )
  rs
  expect_equal( nrow( rs ), 0 )

  rs = result_table( umt, return= "sc_units" )
  rs
  expect_equal( nrow(rs), 8 )

  rs = result_table( umt, return= "agg_co_units" ) %>%
    dplyr::filter( Z == 0 )
  rs
  expect_true( all(rs$id %in% as.character( dat$id ) ) )


  # And if we have no matches left?

  expect_warning(
    mtch <- get_cal_matches(dat,
                            form = Z ~ X1 + X2,
                            metric = "maximum",
                            scaling = c(1/0.2, 1/0.2),
                            caliper = 0.05,
                            rad_method = "fixed",
                            est_method = "scm")
  )
  expect_output(
    summary( mtch ),
    "0 treated units matched to 0 of 4 control units"
  )

  tt = CSM:::unmatched_units(mtch)
  expect_equal( nrow(tt), 5 )


})


test_that("result_table schemas are consistent across return modes", {

  set.seed(4044440)
  dat <- gen_one_toy(nt = 5)

  mtch <- get_cal_matches(dat,
                          treatment = "Z",
                          metric = "maximum",
                          scaling = c(1/0.2, 1/0.2),
                          caliper = 0.25,
                          rad_method = "adaptive",
                          est_method = "scm")

  # all table should contain outcome columns from matches
  rt_all <- result_table(mtch, return = "all")
  expect_true(all(c("id","subclass","weights","Z") %in% names(rt_all)))
  expect_true("Y" %in% names(rt_all))   # if gen_one_toy includes Y; if not, adjust to your outcome column

  # sc_units often drops Y unless explicitly requested
  rt_sc <- result_table(mtch, return = "sc_units")
  expect_true(all(c("id","subclass","weights","Z") %in% names(rt_sc)))

  # If your design is: sc_units doesn't include outcome by default:
  expect_false("Y" %in% names(rt_sc))

  # But should include outcome if asked
  rt_sc_y <- result_table(mtch, return = "sc_units", outcome = "Y")
  expect_true("Y" %in% names(rt_sc_y))
})
