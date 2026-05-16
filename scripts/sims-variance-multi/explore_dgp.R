# scripts/sims-variance-multi/explore_dgp.R
#
# DGP sanity check: generate one large dataset per design cell and print
# variance decompositions so we can verify that the simulation scenarios
# are sensibly calibrated.
#
# Cells explored: error_type × sigma1_extra × overlap_label  (5×2×2 = 20)
# k_match is ignored here — it only affects the matching step, not the DGP.
#
# Sample sizes: nt = 2000, nc = 10000  (large enough that sampling noise
# is negligible and population quantities are well estimated).

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(knitr)
})

source(here::here("scripts/sims-variance-multi/0_utils.R"))

# Re-attach dplyr verbs in case sourcing 0_utils.R caused masking
library(dplyr)

NT <- 2000
NC <- 10000
SEED <- 42

# ── Unique DGP cells ─────────────────────────────────────────────────────────

dgp_grid <- DESIGN_GRID %>%
  distinct(overlap_label, error_type, sigma1_extra, common_label, tx_type) %>%
  arrange(error_type, common_label, overlap_label)

cat("Exploring", nrow(dgp_grid), "DGP cells  (nt =", NT, ", nc =", NC, ")\n\n")


# ── Helper: summarise one dataset ────────────────────────────────────────────

summarise_dgp <- function(df) {
  tx  <- df %>% filter(Z == 1)
  ctl <- df %>% filter(Z == 0)

  # noise decomposition: Var(Y) = Var(signal) + E[sigma^2]
  var_Y0_signal <- var(tx$Y0_denoised)
  var_Y1_signal <- var(tx$Y1_denoised)
  cate_i  <- tx$Y1_denoised - tx$Y0_denoised  # individual treatment effects (signal part)
  tau_i         <- tx$Y1 - tx$Y0  # individual treatment effects

  tibble(
    # ── Treated units ────────────────────────────────────────────────────
    nt              = nrow(tx),
    nc              = nrow(ctl),
    # noise SDs
    sig0_tx_mean    = mean(tx$sigma0),
    sig0_tx_sd      = sd(tx$sigma0),
    sig1_tx_mean    = mean(tx$sigma1),
    sig1_tx_sd      = sd(tx$sigma1),
    # noise variances  E[sigma^2]
    E_sig0sq_tx     = mean(tx$sigma0^2),
    E_sig1sq_tx     = mean(tx$sigma1^2),
    # signal variances
    var_Y0sig_tx    = var_Y0_signal,
    var_Y1sig_tx    = var_Y1_signal,
    # total outcome variances
    var_Y0_tx       = var(tx$Y0),
    var_Y1_tx       = var(tx$Y1),
    # treatment effect distribution
    tau_mean        = mean(tau_i),
    tau_sd          = sd(tau_i),
    cate_mean = mean(cate_i),
    cate_sd   = sd(cate_i),
    # ── Control units ────────────────────────────────────────────────────
    sig0_ctl_mean   = mean(ctl$sigma0),
    sig0_ctl_sd     = sd(ctl$sigma0),
    E_sig0sq_ctl    = mean(ctl$sigma0^2),
    var_Y0_ctl      = var(ctl$Y0),
    # ── SNR (signal / noise, treated side) ───────────────────────────────
    snr_Y0_tx       = sqrt(var_Y0_signal) / mean(tx$sigma0),
    snr_tau         = mean(tau_i)         / mean(tx$sigma0)
  )
}


# ── Run over all cells ────────────────────────────────────────────────────────

set.seed(SEED)

results <- pmap_dfr(dgp_grid, function(overlap_label, error_type,
                                        sigma1_extra, common_label) {
  df <- make_df_multi(
    nc           = NC,
    nt           = NT,
    prop_nc_unif = PROP_NC_UNIF[as.character(overlap_label)],
    error_type   = as.character(error_type),
    sigma1_extra = sigma1_extra
  )

  summarise_dgp(df) %>%
    mutate(
      overlap_label = as.character(overlap_label),
      error_type    = as.character(error_type),
      common_label  = as.character(common_label),
      .before       = everything()
    )
})


# ── Print summaries ───────────────────────────────────────────────────────────

cat("══════════════════════════════════════════════════════════════════════\n")
cat("1. NOISE STANDARD DEVIATIONS  (treated units)\n")
cat("══════════════════════════════════════════════════════════════════════\n")
results %>%
  dplyr::select(error_type, common_label, overlap_label,
         sig0_tx_mean, sig0_tx_sd, sig1_tx_mean, sig1_tx_sd) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  kable(col.names = c("error", "common", "overlap",
                       "σ₀ mean", "σ₀ sd", "σ₁ mean", "σ₁ sd"),
        format = "simple") %>%
  print()

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("2. NOISE VARIANCES  E[σ²]  vs  SIGNAL VARIANCE  Var(signal)\n")
cat("   (treated side; Var(Y) ≈ Var(signal) + E[σ²])\n")
cat("══════════════════════════════════════════════════════════════════════\n")
results %>%
  dplyr::select(error_type, common_label, overlap_label,
         E_sig0sq_tx, var_Y0sig_tx, var_Y0_tx,
         E_sig1sq_tx, var_Y1sig_tx, var_Y1_tx) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  kable(col.names = c("error", "common", "overlap",
                       "E[σ₀²]", "Var(f₀)", "Var(Y0)",
                       "E[σ₁²]", "Var(f₁)", "Var(Y1)"),
        format = "simple") %>%
  print()

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("3. TREATMENT EFFECT DISTRIBUTION  (true τᵢ = Y1_signal - Y0_signal)\n")
cat("══════════════════════════════════════════════════════════════════════\n")
results %>%
  dplyr::select(error_type, common_label, overlap_label,
         tau_mean, tau_sd, cate_mean, cate_sd ) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  kable(col.names = c("error", "common", "overlap",
                       "mean(τᵢ)", "sd(τᵢ)", "mean(cate)", "sd(cate)"),
        format = "simple") %>%
  print()

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("4. CONTROL-SIDE NOISE  (σ₀ for controls varies by overlap)\n")
cat("══════════════════════════════════════════════════════════════════════\n")
results %>%
  dplyr::select(error_type, common_label, overlap_label,
         sig0_ctl_mean, sig0_ctl_sd, E_sig0sq_ctl, var_Y0_ctl) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  kable(col.names = c("error", "common", "overlap",
                       "σ₀ mean", "σ₀ sd", "E[σ₀²]", "Var(Y0_ctl)"),
        format = "simple") %>%
  print()

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("5. SIGNAL-TO-NOISE RATIOS  (treated side)\n")
cat("   snr_Y0 = sd(f₀) / mean(σ₀);   snr_tau = mean(τ) / mean(σ₀)\n")
cat("══════════════════════════════════════════════════════════════════════\n")
results %>%
  dplyr::select(error_type, common_label, overlap_label,
         snr_Y0_tx, snr_tau) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  kable(col.names = c("error", "common", "overlap",
                       "SNR(Y0)", "SNR(τ)"),
        format = "simple") %>%
  print()

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("6. WIDE SUMMARY  (all key quantities, for copy-paste)\n")
cat("══════════════════════════════════════════════════════════════════════\n")
results %>%
  dplyr::select(error_type, common_label, overlap_label,
         sig0_tx_mean, sig1_tx_mean,
         E_sig0sq_tx, E_sig1sq_tx,
         var_Y0_tx, var_Y1_tx,
         tau_mean, tau_sd,
         snr_tau) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  print(n = Inf)

# Return invisibly in case called interactively
invisible(results)
