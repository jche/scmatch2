

test_that("aggregation works", {

  mtch = list(
    data.frame( id = c("T1", "C1", "C2"),
                Z = c(1,0,0),
                Y = c(5,2,3),
                X1 = c(1,2,3),
                X2 = c(4,5,6),
                weights = c(1, 0.6, 0.4),
                subclass = c("T1", "T1", "T1") ),
    data.frame( id = c("T2", "C1"),
                Z = c(1,0),
                Y = c(7,4),
                X1 = c(2,3),
                X2 = c(5,7),
                weights = c(1,1),
                subclass = c("T2", "T2") ),
    data.frame( id = c("T3", "C4", "C5", "C6"),
                Z = c(1,0,0,0),
                Y = c(6,2,3,5),
                X1 = c(3,4,5,6),
                X2 = c(7,8,9,10),
                weights = c(1, 0.5, 0.3, 0.2),
                subclass = c("T3", "T3", "T3", "T3") )
    )

  aa <- agg_avg_units( mtch,
                 covariates = c("X1","X2"),
                 treatment = "Z",
                 outcome = "Y" )
  aa
  expect_equal( aa$X1[[5]], mean( mtch[[3]]$X1[ -1 ] ) )

  co <- agg_co_units( mtch,
                      covariates = c("X1","X2"),
                      treatment = "Z",
                      outcome = "Y" )
  co
  expect_equal( co$weights[[1]], 1.6 )

  cc <- agg_sc_units( mtch,
                      covariates = c("X1","X2"),
                      treatment = "Z",
                      outcome = "Y" )
  cc
  X2 = weighted.mean( mtch[[3]]$X2[-1], mtch[[3]]$weights[-1] )
  expect_equal( cc$X2[[5]], X2 )

})
