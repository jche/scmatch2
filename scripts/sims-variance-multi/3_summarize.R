# scripts/sims-variance-multi/3_summarize.R
#
# Compute coverage (with MCSEs via simhelpers) and produce the paper coverage plot.
#
# Plot specification (from paper writeup):
#   x-axis : overlap level
#   y-axis : coverage
#   colour : variance estimator (homo, het, alt_common, alt_tt)
#   facets : error_type (rows) × common_label (columns)
#
# Usage:
#   Rscript 3_summarize.R [output_name]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(here)
  library(simhelpers)    # calc_coverage() with MCSEs
})

args        <- commandArgs(trailingOnly = TRUE)
output_name <- if (length(args) >= 1) args[[1]] else "sims-variance-multi"

source(here::here("scripts/sims-variance-multi/0_utils.R"))
paths <- get_sim_paths(output_name)

# ── Load combined results ────────────────────────────────────────────────────
combined_path <- paths$combined_rds
if (!file.exists(combined_path))
  combined_path <- paths$combined_csv

stopifnot("Run 2_collect.R first" = file.exists(combined_path))

res <- if (grepl("\\.rds$", combined_path)) {
  readRDS(combined_path)
} else {
  read_csv(combined_path, show_col_types = FALSE)
}

cat("Loaded ", nrow(res), " rows from ", basename(combined_path), "\n", sep = "")

# ── Coverage column ──────────────────────────────────────────────────────────
res <- res %>%
  mutate(
    CI_lower = att_est - 2 * SE,
    CI_upper = att_est + 2 * SE,
    covered = (!is.na(CI_lower) & !is.na(CI_upper) & !is.na(SATT)) &
              (CI_lower <= SATT) & (SATT <= CI_upper),
    # ordered factors for clean plots
    overlap_label = factor(overlap_label, levels = OVERLAP_LABELS),
    error_type    = factor(error_type,
                           levels = c("homo", "het"),
                           labels = c("Homoskedastic", "Heteroskedastic")),
    common_label  = factor(common_label,
                           levels = c("common", "no_common"),
                           labels = c("Common variance\n(σ₁ = σ₀)",
                                      "No common variance\n(σ₁ = σ₀ + 2)")),
    method = factor(method,
                              levels = c("homo", "het", "ttmatch" ),
                              labels = c("homo", "het", "ttmatch" ) )
  )


# ── Compute coverage + MCSE by (estimator × design cell) ────────────────────
coverage_tbl <- res %>%
  group_by(method, overlap_label, error_type, common_label) %>%
  reframe(
    simhelpers::calc_coverage(
      pick(everything()),
      lower_bound = CI_lower,
      upper_bound = CI_upper,
      true_param  = SATT
    )
  ) %>%
  rename(coverage = coverage, mcse = coverage_mcse)

cat("\nCoverage summary (first rows):\n")
print(head(coverage_tbl, 20))

# ── Save summary CSV ─────────────────────────────────────────────────────────
out_dir <- here::here("tables", output_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(coverage_tbl, file.path(out_dir, "coverage_table.csv"))
cat("Saved coverage table: ", file.path(out_dir, "coverage_table.csv"), "\n")

# ── Coverage plot ─────────────────────────────────────────────────────────────
overlap_x_labels <- c(
  "very_low"  = "Very Low",
  "low"       = "Low",
  "mid"       = "Medium",
  "high"      = "High",
  "very_high" = "Very High"
)

estimator_colours <- c(
  "homo"    = "#E41A1C",
  "het"     = "#377EB8",
  "alt-c"   = "#4DAF4A",
  "alt-tt"  = "#984EA3"
)

p_coverage <- ggplot(
  coverage_tbl,
  aes(x = overlap_label, y = coverage,
      colour = method, group = method)
) +
  # ±2 MCSE error bars
  geom_linerange(
    aes(ymin = coverage - 2 * mcse,
        ymax = coverage + 2 * mcse),
    linewidth = 0.5, alpha = 0.6,
    position  = position_dodge(width = 0.3)
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  facet_grid(error_type ~ common_label) +
  scale_x_discrete(labels = overlap_x_labels) +
  scale_y_continuous(
    limits = c(0.7, 1.0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_manual(values = estimator_colours, name = "Estimator") +
  labs(
    x     = "Degree of overlap",
    y     = "Coverage of τ_SATT",
    title = "Coverage of 95% CIs across overlap, error structure, and common-variance scenarios"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.text      = element_text(size = 9)
  )

fig_dir <- here::here("figures", output_name)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
fig_path <- file.path(fig_dir, "coverage_by_overlap.pdf")
ggsave(fig_path, p_coverage, width = 8, height = 6)
cat("Saved figure: ", fig_path, "\n")

print(p_coverage)
