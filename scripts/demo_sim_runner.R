# Demo the sim_runner method

library(tidyverse)
library(mvtnorm)
library(optweight)
library(dbarts)
library(here)

devtools::load_all()
source(here::here("scripts/lib/wrappers.R"))
source(here::here("scripts/lib/sim_runner.R"))

set.seed(42)

# ------------------------------------------------------------------------------
# 1. Data generation
# ------------------------------------------------------------------------------
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


cat(sprintf("Dataset: n=%d  (nt=%d, nc=%d)\n", nrow(df), sum(df$Z), sum(!df$Z)))

true_ATT <- df %>%
  filter(Z == 1) %>%
  summarize(att = mean(Y1 - Y0)) %>%
  pull(att)

cat(sprintf("True ATT: %.4f\n\n", true_ATT))

# ------------------------------------------------------------------------------
# 2. Scaling + feasibility filter
# ------------------------------------------------------------------------------

covs <- parse_form(form)$covs
dist_scaling <- df %>%
  summarize(across(all_of(covs), function(x) nbins / (max(x) - min(x))))

preds_csm <- get_cal_matches(data = df, treatment = "Z",
                             metric = "maximum", scaling = dist_scaling,
                             rad_method = "fixed", est_method = "average", k = 25)
preds_cem <- get_cem_matches(data = df, num_bins = nbins, est_method = "average",
                             return = "all")

ninf     <- length(attr(preds_csm, "unmatched_units"))
ninf_cem <- sum(df$Z) - length(attr(preds_cem, "feasible_units"))
df       <- df %>% filter(!id %in% attr(preds_csm, "unmatched_units"))

cat(sprintf("Infeasible units (CSM): %d,  (CEM): %d\n", ninf, ninf_cem))
cat(sprintf("Analysis n after filter: %d  (nt=%d)\n\n", nrow(df), sum(df$Z)))



# ------------------------------------------------------------------------------
# 3. Run all fast methods  (pass "all" for slow ones too)
# ------------------------------------------------------------------------------


results <- run_all_methods(
  df               = df,
  form             = form,
  dist_scaling     = dist_scaling,
  nbins            = nbins,
  selected_methods = "fast"
) %>%
  mutate(
    true_ATT = true_ATT,
    bias     = ATT_est - true_ATT
  )

cat("\n--- Results ---\n")
cat(sprintf("True ATT: %.4f\n\n", true_ATT))
print( results %>%
        dplyr::select(method, ATT_est, bias, secs), n = Inf)




