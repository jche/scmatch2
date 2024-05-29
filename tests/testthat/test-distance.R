

test_that("Making the {{ }} call work", {
  df <- data.frame(X1=1:4,
                   X2=4:1,
                   Z = c(0, 1, 0, 1))
  dm <- CSM:::gen_dm(df,
               covs=starts_with("X"),
               treatment=Z)
  expect_true(is.matrix(dm))
})

test_that("Error when covariates or treatment variable specified in covs or treatment are not in the dataframe",{
  df <- data.frame(X1 = 1:5, Z=rep(0:1, length.out=5))
  expect_error(gen_dm(df, covs=c("X1","X2") , treatment = "Z"))
  expect_error(gen_dm(df, covs=c("X1") , treatment = "A"))
})

test_that("coerce_covs coerces T/F into 1/2, and turn factor types into numerics",{
  covs <- data.frame(X1=as.factor(1:4),
                      X2=c(T,F,F,T))
  expected_coerced_covs <-
    data.frame(X1=1:4,
               X2=c(2,1,1,2))

  actual_coerced_covs <- coerce_covs(covs)

  expect_equal(expected_coerced_covs,
               actual_coerced_covs)
})

test_that("scale_covs correctly scales the input data",{
  # 1-d case
  covs <- data.frame(X1=1:4)
  scaling <- 2

  expected <-
    data.frame(X1=(1:4)*2) %>%
    as.matrix()

  actual <- scale_covs(covs, scaling)

  dimnames(actual) <- NULL
  dimnames(expected) <- NULL

  expect_equal(expected,
               actual)
  expect_true( is.matrix(actual) )

  # 2-d case
  covs <- data.frame(X1=1:4,
                     X2=c(2,1,1,2))
  scaling <- c(2, 1/2)

  expected <-
    data.frame(X1=(1:4)*2,
              X2=c(2,1,1,2)/2) %>%
    as.matrix()

  actual <- scale_covs(covs, scaling)

  dimnames(actual) <- NULL
  dimnames(expected) <- NULL
  expect_true( is.matrix(actual) )

  expect_equal(expected,
               actual)

})

test_that("Manual computation should match the gen_dm output",{
  df <- data.frame(X1=1:4,
                   X2=4:1,
                   Z = c(0, 1, 0, 1))
  covs <- c("X1", "X2")
  treatment <- "Z"
  expected_dm <-
    matrix(c(sqrt(2), sqrt(18), sqrt(2),  sqrt(2)), nrow = 2)

  actual_dm <- CSM:::gen_dm(df, covs, treatment,
                      scaling=1,
                      metric = "euclidean")
  # only check
  dimnames(actual_dm) <- NULL
  dimnames(expected_dm) <- NULL
  expect_equal(as.matrix(actual_dm), expected_dm, ignore_attr = FALSE)
})
