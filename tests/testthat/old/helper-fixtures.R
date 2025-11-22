# tests/testthat/helper-fixtures.R
make_mock_matches <- function() {
  set.seed(123)

  # Small synthetic dataset (use your existing toy generator)
  df <- gen_one_toy(k = 2, nc = 40, nt = 10, f0_sd = 0.1)

  # Simple scaling for CSM (keep tiny to be fast)
  nbins <- 5
  scaling <- df |>
    dplyr::summarize(dplyr::across(dplyr::starts_with("X"),
                                   ~ nbins / (max(.) - min(.))
    ))

  # Produce a matches table/data.frame the SE functions accept
  # (treats and controls with subclass + weights)
  matches <- CSM::get_cal_matches(
    df = df,
    metric = "maximum",
    scaling = scaling,
    rad_method = "fixed",
    est_method = "average",
    k = 10
  )

  # If your SE funcs expect a long table (subclass, Z, Y, weights, id),
  # convert the CSM object; adjust this to your actual helper if different.
  # If you already have a helper, use it. Otherwise:
  # out <- attr(matches, "match_table")
  out <- full_unit_table(matches, nonzero_weight_only = TRUE )
  if (is.null(out)) {
    # Fallback: build from df + matches attributes
    # (Replace this with the correct converter you have in your package)
    stop("Please replace with your converter from CSM matches to long table.")
  }
  out
}
