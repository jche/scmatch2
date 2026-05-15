# scripts/ferman-analysis/04-ferman-caliper-sensitivity.R
#
# Caliper misspecification sensitivity analysis for the Ferman (Brazil) data.
#
# Motivation (reviewer request): show how SATT and FSATT estimates change when
# the analyst makes different covariate-caliper choices, to illustrate that
# CSM results are robust (or characterise their sensitivity) to user choices
# of pi (the scaling / caliper parameters).
#
# Design:
#   y2007, y2008, y2009 calipers each varied at {0.5x, 1x, 2x} of the
#   base value (0.2) → 3^3 = 27 combinations.
#   is_sao_paolo caliper at {1/1000 (exact/current), 2 (loose)} → 2 levels.
#   Total: 27 × 2 = 54 specifications.
#
# Output:
#   - figures/ferman_caliper_caterpillar.pdf  (caterpillar plot)
#   - tables/ferman_caliper_sensitivity.csv   (full results table)

suppressPackageStartupMessages({
  library(CSM)
  library(tidyverse)
  library(here)
  library(patchwork)
})

# ── Load data & base specification ──────────────────────────────────────────
# 02-core-ferman-analysis.R sets:
#   ferman_for_analysis, match_covs, c (global caliper = 0.35),
#   covariate_caliper = c(0.2, 0.2, 0.2, 1/1000), scaling
source(here::here("scripts/ferman-analysis/02-core-ferman-analysis.R"))

BASE_CAL <- covariate_caliper   # c(0.2, 0.2, 0.2, 0.001)
GLOBAL_C <- c                   # 0.35  (held fixed throughout)

# ── Specification grid ───────────────────────────────────────────────────────
year_mults <- c(0.5, 1, 2)                # multipliers for y2007/08/09
sp_cals    <- c(1 / 1000, 2)             # exact match vs. effectively ignored

spec_grid <- expand.grid(
  mult_y2007 = year_mults,
  mult_y2008 = year_mults,
  mult_y2009 = year_mults,
  cal_sp     = sp_cals,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Human-readable spec label: "07:½ 08:1 09:2 SP:E"
mult_label <- function(m) c("0.5" = "½", "1" = "1", "2" = "2")[as.character(m)]
sp_label   <- function(x) ifelse(x < 0.1, "E", "L")   # E=exact, L=loose

spec_grid <- spec_grid %>%
  mutate(
    spec_id = row_number(),
    label   = paste0(
      mult_label(mult_y2007), "-",
      mult_label(mult_y2008), "-",
      mult_label(mult_y2009), "-",
      sp_label(cal_sp)
    ),
    is_base = (mult_y2007 == 1 & mult_y2008 == 1 & mult_y2009 == 1 & cal_sp < 0.1)
  )

cat("Running", nrow(spec_grid), "caliper specifications...\n")


# ── Run one specification ────────────────────────────────────────────────────
run_one_spec <- function(mult_y2007, mult_y2008, mult_y2009, cal_sp,
                         spec_id, label, is_base) {

  cal_vec <- c(
    BASE_CAL[1] * mult_y2007,
    BASE_CAL[2] * mult_y2008,
    BASE_CAL[3] * mult_y2009,
    cal_sp
  )
  sc <- 1 / cal_vec

  csm_fit <- tryCatch(
    ferman_for_analysis %>%
      get_cal_matches(
        covs       = match_covs,
        treatment  = "Z",
        caliper    = GLOBAL_C,
        metric     = "maximum",
        rad_method = "adaptive",
        scaling    = sc,
        est_method = "scm"
      ),
    error = function(e) {
      message("spec ", spec_id, " (", label, "): matching failed — ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(csm_fit)) return(NULL)

  extract_est <- function(feasible_only, estimand) {
    tryCatch({
      res <- estimate_ATT(csm_fit, outcome = "Y", feasible_only = feasible_only)
      tibble(
        estimand  = estimand,
        ATT       = res$ATT,
        SE        = res$SE,
        N_T       = res$N_T,
        ESS_C     = res$ESS_C
      )
    }, error = function(e) {
      message("spec ", spec_id, " ", estimand, " failed: ", conditionMessage(e))
      NULL
    })
  }

  bind_rows(
    extract_est(feasible_only = FALSE, estimand = "SATT"),
    extract_est(feasible_only = TRUE,  estimand = "FSATT")
  ) %>%
    mutate(
      spec_id    = spec_id,
      label      = label,
      is_base    = is_base,
      mult_y2007 = mult_y2007,
      mult_y2008 = mult_y2008,
      mult_y2009 = mult_y2009,
      cal_sp     = cal_sp,
      sp_exact   = cal_sp < 0.1
    )
}

# ── Run all specs ────────────────────────────────────────────────────────────
results_list <- pmap(spec_grid, run_one_spec)
results <- bind_rows(results_list)

cat("Completed", nrow(results) / 2, "specs with both SATT & FSATT\n")

# Save raw results
tab_dir <- here::here("tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(results, file.path(tab_dir, "ferman_caliper_sensitivity.csv"))
cat("Saved results table\n")


# ── Caterpillar plot ─────────────────────────────────────────────────────────
# Sort specs by ATT within each estimand; base spec highlighted.

plot_df <- results %>%
  mutate(
    ci_lo = ATT - 2 * SE,
    ci_hi = ATT + 2 * SE,
    label = factor(label)
  ) %>%
  mutate(label = fct_reorder(label, ATT))

base_line <- results %>%
  filter(is_base) %>%
  select(estimand, ATT) %>%
  rename(base_ATT = ATT)

plot_df <- left_join(plot_df, base_line, by = "estimand")

p_cat <- ggplot(plot_df,
                aes(x = ATT, y = label,
                    colour = is_base, size = is_base)) +
  geom_vline(aes(xintercept = base_ATT),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0, linewidth = 0.4, alpha = 0.6) +
  geom_point() +
  scale_colour_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey40"),
                      guide = "none") +
  scale_size_manual(values = c("TRUE" = 2.5, "FALSE" = 1.2),
                    guide = "none") +
  facet_wrap(~ estimand) +
  labs(
    x        = "ATT estimate (±2 SE)",
    y        = "Caliper spec  [y07-y08-y09-SP]",
    title    = "Sensitivity to covariate caliper choices (Ferman data)",
    subtitle = "Red = base specification.  Dashed line = base ATT.  SP: E = exact match, L = loose."
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text         = element_text(size = 9),
    axis.text.y        = element_text(size = 7),
    panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90")
  )

print(p_cat)

fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
fig_path <- file.path(fig_dir, "ferman_caliper_caterpillar.pdf")
ggsave(fig_path, p_cat, width = 6, height = 8)
cat("Saved figure:", fig_path, "\n")



plot_df

dfW <- plot_df %>%
  dplyr::select( label, ATT, estimand ) %>%
  pivot_wider( names_from = estimand, values_from = ATT )

plot2 <- ggplot( dfW, aes( SATT, FSATT ) ) +
  geom_point() +
  geom_abline( slope=1, intercept=0, linetype="dashed", colour="grey50" ) +
  labs( x = "SATT", y = "FSATT" ) +
  theme_minimal() +
  coord_fixed()

fig_path <- file.path(fig_dir, "ferman_caliper_FSATTvSATT.pdf")
ggsave(fig_path, plot2, width = 5, height = 5)
cat("Saved figure:", fig_path, "\n")


