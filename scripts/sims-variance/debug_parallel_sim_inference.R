source(here::here("scripts/sims-variance/0_sim_inference_utils.R"))

overlap_label <- "very_low"
error_label <- "homoskedastic"

prop_nc_unif_values <- c(
  very_low = 2/3,
  low = 1/2,
  mid = 1/3,
  high = 1/5,
  very_high = 1/10
)
prop_nc_unif <- prop_nc_unif_values[[overlap_label]]

set.seed(123)

cat("\n====================\n")
cat("DEBUG: parameters\n")
cat("====================\n")
print(list(
  overlap_label = overlap_label,
  error_label = error_label,
  prop_nc_unif = prop_nc_unif,
  nc = 500,
  nt = 100,
  f0_sd = 0.5,
  ctr_dist = 0.5,
  nbins = 6
))

cat("\n====================\n")
cat("DEBUG: build df (gen_df_adv)\n")
cat("====================\n")

# replicate your intended DGP
nc <- 500; nt <- 100; f0_sd <- 0.5
df <- gen_df_adv(
  nc = nc, nt = nt,
  f0_sd = f0_sd,
  ctr_dist = 0.5,
  prop_nc_unif = prop_nc_unif,
  tx_effect_fun = function(X1, X2) { 3*X1 + 3*X2 },
  f0_fun = function(x, y) {
    matrix(c(x,y), ncol=2) %>%
      mvtnorm::dmvnorm(
        mean = c(0.5, 0.5),
        sigma = matrix(c(1, 0.8, 0.8, 1), 2)
      ) * 20
  }
)

cat("DEBUG: df nrow =", nrow(df), "\n")
cat("DEBUG: df columns:\n")
print(names(df))

cat("\nDEBUG: df summary (key cols)\n")
print(df %>% dplyr::summarise(
  n = dplyr::n(),
  n_treated = sum(Z == 1),
  n_control = sum(Z == 0),
  Y_mean = mean(Y),
  Y_sd = sd(Y),
  X1_rng = paste0("[", round(min(X1),3), ", ", round(max(X1),3), "]"),
  X2_rng = paste0("[", round(min(X2),3), ", ", round(max(X2),3), "]")
))

cat("\nDEBUG: SATT check (true sample ATT)\n")
true_satt <- df %>%
  dplyr::filter(Z == 1) %>%
  dplyr::summarise(att = mean(Y1 - Y0)) %>%
  dplyr::pull(att)
print(true_satt)

cat("\n====================\n")
cat("DEBUG: scaling used by matching\n")
cat("====================\n")

nbins <- 6
dist_scaling <- df %>%
  dplyr::summarise(dplyr::across(
    dplyr::starts_with("X"),
    function(x) (2*nbins) / (max(x) - min(x))
  )) %>%
  as.numeric()

cat("DEBUG: dist_scaling =\n")
print(dist_scaling)

cat("\n====================\n")
cat("DEBUG: run matching (get_cal_matches)\n")
cat("====================\n")

mtch <- tryCatch(
  get_cal_matches(
    data = df,
    form = Z ~ X1 + X2,
    metric = "maximum",
    scaling = dist_scaling,
    # IMPORTANT: use adaptive caliper like your paper text
    rad_method = "adaptive",
    caliper = 0.2,
    # choose estimator consistent with your pipeline
    est_method = "scm",
    k = 5,
    warn = FALSE
  ),
  error = function(e) e
)

if (inherits(mtch, "error")) {
  cat("ERROR: matching failed:\n")
  print(mtch$message)
  stop(mtch)
}

cat("DEBUG: matching succeeded. class(mtch) =\n")
print(class(mtch))
cat("DEBUG: params(mtch)$treatment =\n")
print(params(mtch)$treatment)

cat("\n====================\n")
cat("DEBUG: result tables schema\n")
cat("====================\n")

rt_all <- result_table(mtch, return = "all")
cat("DEBUG result_table(mtch,'all') names:\n")
print(names(rt_all))
cat("DEBUG nrow(all) =", nrow(rt_all), "\n")
cat("DEBUG unique Z in all:\n")
print(unique(rt_all$Z))
cat("DEBUG has weights/subclass/id/Y?\n")
print(list(
  has_weights = "weights" %in% names(rt_all),
  has_subclass = "subclass" %in% names(rt_all),
  has_id = "id" %in% names(rt_all),
  has_Y = "Y" %in% names(rt_all),
  has_Z = "Z" %in% names(rt_all)
))

rt_sc <- tryCatch(result_table(mtch, return = "sc_units"), error = function(e) e)
if (!inherits(rt_sc, "error")) {
  cat("\nDEBUG result_table(mtch,'sc_units') names:\n")
  print(names(rt_sc))
  cat("DEBUG nrow(sc_units) =", nrow(rt_sc), "\n")
  cat("DEBUG has Y in sc_units? ", "Y" %in% names(rt_sc), "\n")
} else {
  cat("\nERROR: result_table(mtch,'sc_units') failed:\n")
  print(rt_sc$message)
}

cat("\n====================\n")
cat("DEBUG: try get_ATT_estimate\n")
cat("====================\n")

att_try <- tryCatch(
  get_ATT_estimate(
    mtch,
    treatment = "Z",
    outcome = "Y",
    variance_method = "pooled"
  ),
  error = function(e) e
)

if (inherits(att_try, "error")) {
  cat("\nERROR: get_ATT_estimate failed\n")
  cat("Message:", att_try$message, "\n")
  cat("\nTRACEBACK:\n")
  tb <- sys.calls()
  print(tb)
} else {
  cat("\nDEBUG: get_ATT_estimate succeeded:\n")
  print(att_try)
}

cat("\n====================\n")
cat("DEBUG: isolate ATT point est vs variance\n")
cat("====================\n")

att_point <- tryCatch(get_att_point_est(rt_all, treatment="Z", outcome="Y"), error=function(e) e)
cat("DEBUG: get_att_point_est(rt_all) =>\n")
print(att_point)

var_only <- tryCatch(get_total_variance(mtch, treatment="Z", outcome="Y", variance_method="pooled"), error=function(e) e)
cat("\nDEBUG: get_total_variance(mtch) =>\n")
print(var_only)

cat("\n====================\n")
cat("DONE\n")
cat("====================\n")
