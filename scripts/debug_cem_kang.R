## debug_cem_kang.R
##
## Diagnose why CEM has high variance on the Kang DGP.
## Hypothesis: the heavily skewed / transformed covariates (X1–X4) create
## strata with extreme tx/co imbalance or empty cells, driving up variance.

devtools::load_all()
suppressPackageStartupMessages({
  library(tidyverse)
  library(mvtnorm)
  library(MatchIt)
})
source(here::here("scripts/lib/wrappers.R"))
source(here::here("scripts/lib/sim_runner.R"))
source(here::here("R/sim_data.R"))

set.seed(4291)

# ── 1. Generate one Kang dataset ─────────────────────────────────────────────
n    <- 1000
nbins <- 5
df   <- gen_df_kang(n = n)
covs <- df |> select(starts_with("X")) |> colnames()
form <- as.formula(paste0("Z ~ ", paste(covs, collapse = "+")))

cat(sprintf("N = %d  |  Treated = %d  |  Control = %d\n",
            nrow(df), sum(df$Z), sum(!df$Z)))



# ── 2. Quick look at covariate distributions ──────────────────────────────────
# The Kang X's are heavily skewed transforms of standard normals;
# this is what CEM has to bin.
cat("\nCovariate summaries (note the skew and scale differences):\n")
df |>
  select(Z, all_of(covs)) |>
  pivot_longer(-Z) |>
  group_by(name, Z) |>
  summarise(
    min  = round(min(value), 2),
    p25  = round(quantile(value, 0.25), 2),
    med  = round(median(value), 2),
    p75  = round(quantile(value, 0.75), 2),
    max  = round(max(value), 2),
    .groups = "drop"
  ) |>
  mutate(Z = if_else(Z, "Treated", "Control")) |>
  print(n = Inf)



# ── 3. Run CEM via MatchIt ────────────────────────────────────────────────────
m_out <- matchit(
  form,
  data      = df,
  method    = "cem",
  estimand  = "ATT",
  cutpoints = nbins,
  k2k       = FALSE
)

table(df$Z )
m_data <- match.data(m_out)
m_data
table( m_data$Z, m_data$subclass )
table( m_data$Z )
cat(sprintf("\nAfter CEM: %d treated and %d control units retained (out of %d tx, %d co)\n",
            sum(m_data$Z), sum(!m_data$Z), sum(df$Z), sum(!df$Z)))


# ── 4. Stratum-level tx / co counts ──────────────────────────────────────────
strata_counts <- m_data |>
  count(subclass, Z) |>
  pivot_wider(names_from = Z, values_from = n,
              names_prefix = "n_", values_fill = 0) |>
  rename(n_co = `n_FALSE`, n_tx = `n_TRUE`) |>
  mutate(
    ratio_co_tx = round(n_co / n_tx, 2),
    # weight each tx unit gets for its matched controls
    weight_per_tx = n_co / n_tx
  ) |>
  arrange(desc(n_tx))

cat("\n── Stratum table (sorted by # treated) ──\n")
print(strata_counts, n = Inf)

# Summary of the imbalance
cat("\n── Stratum summary ──\n")
cat(sprintf("Total strata               : %d\n",  nrow(strata_counts)))
cat(sprintf("Strata with 0 controls     : %d\n",  sum(strata_counts$n_co == 0)))
cat(sprintf("Strata with 1 control      : %d\n",  sum(strata_counts$n_co == 1)))
cat(sprintf("Strata with >=10 controls  : %d\n",  sum(strata_counts$n_co >= 10)))
cat(sprintf("Max controls in one stratum: %d\n",  max(strata_counts$n_co)))
cat(sprintf("Max treated in one stratum : %d\n",  max(strata_counts$n_tx)))


# ── 5. Distribution of per-tx-unit effective weights ─────────────────────────
# Each tx unit gets an ATT contribution weighted by the average of its n_co controls.
# A stratum with 1 tx and 50 controls has the same weight as 1 tx with 1 control,
# BUT the residual variance for that tx unit's cell estimate is huge if n_co is small.

# Effective number of controls matched to each tx unit (= n_co in that stratum)
tx_level <- m_data |>
  filter(Z == TRUE) |>
  left_join(strata_counts |> select(subclass, n_co, n_tx), by = "subclass")

cat("\nDistribution of # controls matched to each tx unit:\n")
tx_level |>
  count(n_co) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print(n = Inf)

cat(sprintf("\nMedian controls per tx unit: %.1f\n", median(tx_level$n_co)))
cat(sprintf("Mean   controls per tx unit: %.1f\n", mean(tx_level$n_co)))

# ── 6. Covariate distribution plots ──────────────────────────────────────────
p_covs <- df |>
  select(Z, all_of(covs)) |>
  pivot_longer(-Z) |>
  ggplot(aes(x = value, fill = factor(Z, labels = c("Control", "Treated")))) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.5) +
  facet_wrap(~name, scales = "free") +
  labs(
    title = "Kang DGP: covariate distributions by treatment status",
    subtitle = "Heavy skew / different supports → uneven CEM cells",
    fill = NULL, x = NULL, y = "Count"
  ) +
  theme_classic() +
  theme(legend.position = "top")
print(p_covs)

# ── 7. Stratum size plot ──────────────────────────────────────────────────────
p_strata <- strata_counts |>
  ggplot(aes(x = n_co, y = n_tx)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "CEM strata: # treated vs # controls",
    subtitle = "Points far above the diagonal = tx-heavy strata (high variance contribution)",
    x = "Controls in stratum", y = "Treated in stratum"
  ) +
  theme_classic()
print(p_strata)

# ── 8. Naive ATT estimate from CEM ───────────────────────────────────────────
# Sanity check: reproduce the CEM ATT estimate by hand
cem_att <- m_data |>
  group_by(subclass) |>
  summarise(
    cell_effect = mean(Y[Z]) - mean(Y[!Z]),
    n_tx        = sum(Z),
    .groups     = "drop"
  ) |>
  summarise(ATT = weighted.mean(cell_effect, n_tx)) |>
  pull(ATT)

cat(sprintf("\nHand-computed CEM ATT estimate: %.4f  (true ATT = 0)\n", cem_att))


# ── 9. Run full method suite on the same dataset ─────────────────────────────
# Mirrors the kang block in sims-bias_mse/run_single_iteration.R:
#   - scaling:  nbins / (max - min)  [no 2× factor; kang is not the toy DGP]
#   - true_ATT: 0  [known by construction for the Kang DGP]
#   - feasibility filter via CSM adaptive caliper (k=25 fixed), as in the sim

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("Running full method suite on the Kang dataset\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

dist_scaling <- df %>%
  summarize(across(
    all_of(covs),
    function(x) if (is.numeric(x)) nbins / (max(x) - min(x)) else 1000
  ))

preds_csm <- get_cal_matches(
  data       = df,
  metric     = "maximum",
  scaling    = dist_scaling,
  rad_method = "fixed",
  est_method = "average",
  k          = 25
)
preds_cem <- get_cem_matches(data = df, num_bins = nbins,
                             est_method = "average", return = "all")

ninf     <- length(attr(preds_csm, "unmatched_units"))
ninf_cem <- sum(df$Z) - length(attr(preds_cem, "feasible_units"))
df_filt  <- df %>% filter(!id %in% attr(preds_csm, "unmatched_units"))

cat(sprintf("Infeasible units — CSM: %d  |  CEM: %d\n", ninf, ninf_cem))
cat(sprintf("Analysis dataset after CSM filter: n=%d  (nt=%d, nc=%d)\n\n",
            nrow(df_filt), sum(df_filt$Z), sum(!df_filt$Z)))

true_ATT <- 0   # Kang DGP has ATT = 0 by construction

results <- run_all_methods(
  df               = df_filt,
  form             = form,
  dist_scaling     = dist_scaling,
  nbins            = nbins,
  selected_methods = "fast",
  verbose          = TRUE
) %>%
  mutate(
    true_ATT = true_ATT,
    bias     = ATT_est - true_ATT
  )

cat("\n── Results (true ATT = 0) ──\n")
print(results %>% select(method, ATT_est, bias, secs), n = Inf)
