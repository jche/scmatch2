# Trying same thing with twang ----

if ( FALSE ) {

  # twang GBM propensity score model
  ps_twang <- ps(
    formula = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8,
    data = as.data.frame( lalonde_df ),
    n.trees = 5000,
    interaction.depth = 3,
    shrinkage = 0.01,
    perm.test.iters = 0,
    stop.method = c("es.mean", "ks.max"),
    estimand = "ATE",
    verbose = FALSE
  )

  # pick stopping rule (use es.mean by default)
  stop_rule <- "es.mean"

  # get propensity scores
  lalonde_df$pscore <- twang::get.weights(ps_twang, stop.method = stop_rule)

  # covariate importance (relative influence from GBM)
  var_imp <- summary(ps_twang, stop.method = stop_rule)$rel.inf %>%
    as_tibble() %>%
    rename(var = var, rel_inf = rel.inf) %>%
    arrange(desc(rel_inf))

  var_imp

  # optional: plot balance + influence
  plot(ps_twang, plots = 1, stop.method = stop_rule)  # balance
  plot(ps_twang, plots = 3, stop.method = stop_rule)  # relative influence

}

