# scripts/debug_distances_kang.R
#
# Compare covariate distances from each treated unit to its matched controls
# under CSM (adaptive caliper + SCM weights) vs. CEM, on the Kang DGP.
#
# For CSM: uses get_distance_table(csm, long_table = FALSE) which returns the
#   SCM-weighted, equally-weighted, and closest-neighbour distances per tx unit.
#
# For CEM: calc_distances_for_CEM() (defined below) mirrors get_distance_table's
#   logic: for each treated unit it computes the average covariate distance to
#   every control unit in the same stratum, using the same scaling and metric.

devtools::load_all()
suppressPackageStartupMessages({
  library(tidyverse)
  library(MatchIt)
  library(here)
})
source(here::here("scripts/lib/wrappers.R"))
source(here::here("R/sim_data.R"))

set.seed(4291)

# ── 1. Generate Kang dataset ─────────────────────────────────────────────────
n     <- 1000
nbins <- 5

df   <- gen_df_kang(n = n)
covs <- df |> dplyr::select(starts_with("X")) |> colnames()
form <- as.formula(paste0("Z ~ ", paste(covs, collapse = "+")))

cat(sprintf("N = %d  |  Treated = %d  |  Control = %d\n\n",
            nrow(df), sum(df$Z), sum(!df$Z)))


# ── 2. Scaling (sim-runner convention: nbins / range) ────────────────────────
dist_scaling <- df %>%
  summarize(across(
    all_of(covs),
    function(x) if (is.numeric(x)) nbins / (max(x) - min(x)) else 1000
  ))

cat("Distance scaling (nbins / range per covariate):\n")
print(dist_scaling)
cat("\n")


# ── 3. Run CSM ───────────────────────────────────────────────────────────────
cat("── Running CSM (adaptive, k=2, scm) ──\n")
csm <- get_cal_matches(
  data       = df,
  formula    = form,
  metric     = "maximum",
  scaling    = dist_scaling,
  rad_method = "adaptive",
  k          = 2,
  est_method = "scm",
  warn       = FALSE
)

n_unmatched_csm <- length(attr(csm, "unmatched_units"))
cat(sprintf("CSM: %d treated matched, %d unmatched\n\n",
            sum(df$Z) - n_unmatched_csm, n_unmatched_csm))


# ── 4. Run CEM via MatchIt ───────────────────────────────────────────────────
cat("── Running CEM (nbins =", nbins, ") ──\n")
m_out  <- matchit(
  form,
  data      = df,
  method    = "cem",
  estimand  = "ATT",
  cutpoints = nbins,
  k2k       = FALSE
)
m_data <- match.data(m_out)

cat(sprintf("CEM: %d treated and %d controls retained (of %d tx, %d co)\n\n",
            sum(m_data$Z), sum(!m_data$Z), sum(df$Z), sum(!df$Z)))


# ── 5. CSM distance table ────────────────────────────────────────────────────
# Returns one row per treated unit with columns: id, SCM, average, closest
cat("── CSM distance table (get_distance_table) ──\n")
csm_dists <- get_distance_table(csm, long_table = FALSE)
cat("CSM distances (first 10 rows):\n")
print(head(csm_dists, 10))
cat(sprintf("\nCSM average distance: mean=%.4f  median=%.4f\n\n",
            mean(csm_dists$average, na.rm = TRUE),
            median(csm_dists$average, na.rm = TRUE)))


# ── 6. CEM distance function ────────────────────────────────────────────────
#
# calc_distances_for_CEM(m_data, covs, scaling, metric = "maximum")
#
# For each CEM stratum, computes the (scaled) covariate distance from every
# treated unit to every control unit in the same stratum, then returns the
# average over controls for each treated unit.
#
# Mirrors get_distance_table's use of gen_dm(): gen_dm() returns an
# (n_tx x n_co) distance matrix for a data frame containing both groups.
# Called per-stratum here because CEM can have multiple treated units per cell.
#
# @param m_data   Output of MatchIt::match.data() — must contain columns:
#                 'id', 'Z' (0/1 or T/F), 'subclass', and all covariate columns.
# @param covs     Character vector of covariate names (same as used for matching).
# @param scaling  Named numeric row/vector of per-covariate scaling factors
#                 (output of the summarize(across(...)) call above).
# @param metric   Distance metric passed to gen_dm() (default "maximum").
#
# @return Tibble with one row per matched treated unit:
#   id, subclass, n_co (controls in stratum), avg_dist (mean distance to controls).

calc_distances_for_CEM <- function(m_data,
                                   covs,
                                   scaling,
                                   metric = "maximum") {

  # Only iterate over strata that contain at least one treated unit
  tx_subclasses <- unique(m_data$subclass[as.logical(m_data$Z)])

  result_list <- vector("list", length(tx_subclasses))

  for (i in seq_along(tx_subclasses)) {
    sc       <- tx_subclasses[[i]]
    sub_data <- m_data %>% filter(subclass == sc)
    tx_rows  <- sub_data %>% filter(as.logical(Z))
    co_rows  <- sub_data %>% filter(!as.logical(Z))

    # Skip degenerate strata (shouldn't happen after match.data, but be safe)
    if (nrow(tx_rows) == 0 || nrow(co_rows) == 0) next

    # gen_dm() returns an (n_tx × n_co) scaled distance matrix.
    # Rows are treated units in the order they appear in sub_data;
    # columns are control units — same ordering as co_rows.
    dm <- gen_dm(
      data      = sub_data,
      covs      = covs,
      treatment = "Z",
      scaling   = scaling,
      metric    = metric
    )

    # rowMeans: average distance from each tx unit to all controls in stratum
    avg_dists <- rowMeans(dm)   # length = nrow(tx_rows)

    result_list[[i]] <- tibble(
      id       = tx_rows$id,
      subclass = sc,
      n_co     = nrow(co_rows),
      avg_dist = avg_dists
    )
  }

  bind_rows(result_list)
}


# ── 7. Compute CEM distances ─────────────────────────────────────────────────
cat("── CEM distance table (calc_distances_for_CEM) ──\n")
cem_dists <- calc_distances_for_CEM(m_data, covs = covs,
                                    scaling = dist_scaling,
                                    metric  = "maximum")
cat("CEM distances (first 10 rows):\n")
print(head(cem_dists, 10))
cat(sprintf("\nCEM average distance: mean=%.4f  median=%.4f  (n_tx = %d)\n\n",
            mean(cem_dists$avg_dist, na.rm = TRUE),
            median(cem_dists$avg_dist, na.rm = TRUE),
            nrow(cem_dists)))


# ── 8. Side-by-side comparison ───────────────────────────────────────────────
cat("── Summary comparison ──\n")
cat(sprintf("%-12s  %8s  %8s  %8s\n", "Method", "Mean", "Median", "N_tx"))
cat(sprintf("%-12s  %8.4f  %8.4f  %8d\n", "CSM (avg)",
            mean(csm_dists$average,   na.rm = TRUE),
            median(csm_dists$average, na.rm = TRUE),
            sum(!is.na(csm_dists$average))))
cat(sprintf("%-12s  %8.4f  %8.4f  %8d\n", "CSM (SCM)",
            mean(csm_dists$SCM,   na.rm = TRUE),
            median(csm_dists$SCM, na.rm = TRUE),
            sum(!is.na(csm_dists$SCM))))
cat(sprintf("%-12s  %8.4f  %8.4f  %8d\n", "CEM (avg)",
            mean(cem_dists$avg_dist,   na.rm = TRUE),
            median(cem_dists$avg_dist, na.rm = TRUE),
            nrow(cem_dists)))


# ── 9. Distribution plot ─────────────────────────────────────────────────────
plot_df <- bind_rows(
  csm_dists %>%
    transmute(id, dist = average, method = "CSM (avg)"),
  csm_dists %>%
    transmute(id, dist = SCM,     method = "CSM (SCM)"),
  cem_dists %>%
    transmute(id, dist = avg_dist, method = "CEM (avg)")
)

p <- ggplot(plot_df, aes(x = dist, fill = method, colour = method)) +
  geom_density(alpha = 0.3, linewidth = 0.7) +
  labs(
    title    = "Distribution of tx-to-control distances: CSM vs. CEM (Kang DGP)",
    subtitle = sprintf("Kang n=%d, nbins=%d, scaling = nbins/range, metric = maximum",
                       n, nbins),
    x        = "Average distance (scaled covariates)",
    y        = "Density",
    fill     = "Method", colour = "Method"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(p)
