# tests/testthat/test-sim_runner.R
source(here::here("scripts/lib/wrappers.R"))
source(here::here("scripts/lib/sim_runner.R"))

# ---------------------------------------------------------------------------
# expand_method_groups() — pure-logic tests, no package dependencies
# ---------------------------------------------------------------------------

test_that("ALL_METHODS and SLOW_METHODS are non-empty character vectors", {
  expect_type(ALL_METHODS, "character")
  expect_type(SLOW_METHODS, "character")
  expect_gt(length(ALL_METHODS), 0)
  expect_gt(length(SLOW_METHODS), 0)
  expect_true(all(SLOW_METHODS %in% ALL_METHODS))
})

test_that("expand_method_groups('all') returns ALL_METHODS in order", {
  expect_equal(expand_method_groups("all"), ALL_METHODS)
})

test_that("expand_method_groups('fast') excludes every slow method", {
  result <- expand_method_groups("fast")
  expect_true(length(result) > 0)
  expect_equal(result, setdiff(ALL_METHODS, SLOW_METHODS))
  expect_false(any(SLOW_METHODS %in% result))
})

test_that("expand_method_groups('slow') returns exactly SLOW_METHODS", {
  result <- expand_method_groups("slow")
  expect_setequal(result, SLOW_METHODS)
})

test_that("expand_method_groups handles named group shorthands", {
  expect_equal(expand_method_groups("matching"),
               c("csm_scm", "csm_avg", "cem_scm", "cem_avg", "onenn"))
  expect_equal(expand_method_groups("ps"),   c("ps_lm", "ps_bart"))
  expect_equal(expand_method_groups("or"),   c("or_lm", "or_bart"))
  expect_equal(expand_method_groups("bal"),  c("bal1", "bal2"))
  expect_equal(expand_method_groups("balance"), c("bal1", "bal2"))
  expect_equal(expand_method_groups("dr"),
               c("tmle1", "aipw1", "tmle2", "aipw2"))
  expect_equal(expand_method_groups("doubly_robust"),
               c("tmle1", "aipw1", "tmle2", "aipw2"))
})

test_that("expand_method_groups passes through individual method names", {
  expect_equal(expand_method_groups("diff"), "diff")
  expect_equal(expand_method_groups(c("diff", "or_lm")), c("diff", "or_lm"))
})

test_that("expand_method_groups deduplicates repeated names/groups", {
  expect_equal(expand_method_groups(c("diff", "diff")), "diff")
  # 'fast' + 'diff' should not repeat 'diff'
  result <- expand_method_groups(c("fast", "diff"))
  expect_equal(result, expand_method_groups("fast"))  # diff already in fast
})

test_that("expand_method_groups warns on unknown method and returns empty", {
  expect_warning(
    result <- expand_method_groups("bogus_method"),
    "Unknown"
  )
  expect_length(result, 0)
})

test_that("expand_method_groups drop_slow removes slow methods", {
  result_all_drop  <- expand_method_groups("all",  drop_slow = TRUE)
  result_fast      <- expand_method_groups("fast")
  result_slow_drop <- expand_method_groups("slow", drop_slow = TRUE)

  expect_equal(result_all_drop, result_fast)
  expect_length(result_slow_drop, 0)
})

test_that("expand_method_groups is case-insensitive for group keys", {
  expect_equal(expand_method_groups("ALL"),     expand_method_groups("all"))
  expect_equal(expand_method_groups("Fast"),    expand_method_groups("fast"))
  expect_equal(expand_method_groups("MATCHING"), expand_method_groups("matching"))
})

# ---------------------------------------------------------------------------
# run_all_methods() — integration test with fast methods only, small dataset
# ---------------------------------------------------------------------------

test_that("run_all_methods returns correct tibble structure (fast methods)", {
  skip_if_not_installed("optweight")

  set.seed(42)
  df <- CSM:::gen_one_toy(nc = 80, nt = 25, f0_sd = 0.5)
  form   <- Z ~ X1 + X2
  nbins  <- 6
  covs   <- parse_form(form)$covs

  dist_scaling <- df %>%
    dplyr::summarize(dplyr::across(
      dplyr::all_of(covs),
      function(x) nbins / (max(x) - min(x))
    ))

  # Standard feasibility filter
  preds_csm <- get_cal_matches(
    data = df, treatment = "Z", metric = "maximum",
    scaling = dist_scaling, rad_method = "fixed",
    est_method = "average", k = 25
  )
  df <- df %>% dplyr::filter(!id %in% attr(preds_csm, "unmatched_units"))

  result <- run_all_methods(
    df               = df,
    form             = form,
    dist_scaling     = dist_scaling,
    nbins            = nbins,
    selected_methods = "fast",
    verbose          = FALSE
  )

  # Shape
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("method", "ATT_est", "secs"))
  expect_equal(nrow(result), length(ALL_METHODS))

  # Fast methods have numeric, non-NA estimates
  fast_methods <- expand_method_groups("fast")
  fast_rows    <- result[result$method %in% fast_methods, ]
  expect_true(all(is.numeric(fast_rows$ATT_est)))
  expect_false(any(is.na(fast_rows$ATT_est)))

  # Slow methods have NA estimates (since we only ran "fast")
  slow_rows <- result[result$method %in% SLOW_METHODS, ]
  expect_true(all(is.na(slow_rows$ATT_est)))
})



test_that("run_all_methods accepts method group shorthands via selected_methods", {
  skip_if_not_installed("optweight")

  set.seed(7)
  df <- CSM:::gen_one_toy(nc = 60, nt = 20, f0_sd = 0.5)
  form  <- Z ~ X1 + X2
  nbins <- 6
  covs  <- parse_form(form)$covs

  dist_scaling <- df %>%
    dplyr::summarize(dplyr::across(
      dplyr::all_of(covs),
      function(x) nbins / (max(x) - min(x))
    ))

  preds_csm <- get_cal_matches(
    data = df, treatment = "Z", metric = "maximum",
    scaling = dist_scaling, rad_method = "fixed",
    est_method = "average", k = 25
  )
  df <- df %>% dplyr::filter(!id %in% attr(preds_csm, "unmatched_units"))

  # Run only the matching methods
  result <- run_all_methods(
    df               = df,
    form             = form,
    dist_scaling     = dist_scaling,
    nbins            = nbins,
    selected_methods = "matching",
    verbose          = FALSE
  )

  matching_methods <- expand_method_groups("matching")
  matching_rows    <- result[result$method %in% matching_methods, ]
  other_rows       <- result[!result$method %in% matching_methods, ]

  expect_true(all(is.numeric(matching_rows$ATT_est)))
  expect_false(any(is.na(matching_rows$ATT_est)))
  expect_true(all(is.na(other_rows$ATT_est)))
})


