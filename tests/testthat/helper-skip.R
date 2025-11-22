# tests/testthat/helper-skip.R
run_slow_tests <- function() {
  isTRUE(as.logical(Sys.getenv("RUN_SLOW_TESTS", "FALSE")))
}
