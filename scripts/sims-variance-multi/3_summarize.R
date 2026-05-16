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

fig_dir <- here::here("figures", output_name)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)



args        <- commandArgs(trailingOnly = TRUE)
output_name <- if (length(args) >= 1) args[[1]] else "sims-variance-multi"

source( here::here( "scripts/lib/plot_sim.R" ) )

source(here::here("scripts/sims-variance-multi/0_utils.R"))
paths <- get_sim_paths(output_name)

# ── Load combined results ────────────────────────────────────────────────────
combined_path <- paths$combined_rds
if (!file.exists(combined_path)) {
  combined_path <- paths$combined_csv
}

stopifnot("Run 2_collect.R first" = file.exists(combined_path))

res <- if (grepl("\\.rds$", combined_path)) {
  readRDS(combined_path)
} else {
  read_csv(combined_path, show_col_types = FALSE)
}

cat("Loaded ", nrow(res), " rows from ", basename(combined_path), "\n", sep = "")

# ── Generate Coverage statistics ----

res <- res %>%
  mutate(
    CI_lower = att_est - 2 * SE,
    CI_upper = att_est + 2 * SE,
    covered = (!is.na(CI_lower) & !is.na(CI_upper) & !is.na(SATT)) &
      (CI_lower <= SATT) & (SATT <= CI_upper)
    )

res$tx_type[ is.na(res$tx_type) ] <- "het"

coverage_tbl <- res %>%
  group_by(method, overlap_label, error_type, common_label, tx_type, k_match ) %>%
  summarize(
    my_calc_coverage(
      pick(everything()),
      lower_bound = CI_lower,
      upper_bound = CI_upper,
      true_param  = SATT
    ),
    R = n(), .groups = "drop"
  ) %>%
  rename(mcse = coverage_mcse)

coverage_tbl$R

cat("\nCoverage summary (first rows):\n")
print(head(coverage_tbl, 10))

# ── Save summary CSV ─────────────────────────────────────────────────────────
out_dir <- here::here("tables", output_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(coverage_tbl, file.path(out_dir, "coverage_table.csv"))
cat("Saved coverage table: ", file.path(out_dir, "coverage_table.csv"), "\n")

# ── Coverage plot ─────────────────────────────────────────────────────────────

coverage_tbl <- coverage_tbl %>%
  mutate(     # ordered factors for clean plots
    overlap_label = factor(overlap_label, levels = OVERLAP_LABELS),
    error_type    = factor(error_type,
                           levels = c("homo", "het"),
                           labels = c("Homoskedastic", "Heteroskedastic")),
    common_label  = factor(common_label,
                           levels = c("common", "no_common"),
                           labels = c("Common variance",
                                      "No common variance")),
    k_match = paste0( "k=", k_match )
  )


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
  "ttmatch"   = "#4DAF4A",
  "pop_het" = "#984EA3",
  "pop_homo" = "#FF7F00"
)

coverage_tbl <- filter( coverage_tbl,
                        method != "pop_het" & method != "pop_homo" )
p_coverage <- ggplot(
  filter( coverage_tbl, tx_type == "het" ),
  aes(x = overlap_label, y = coverage,
      colour = method, group = method)
) +
  facet_grid( error_type ~ k_match + common_label ) +
  geom_linerange(
    aes(ymin = coverage - 2 * mcse,
        ymax = coverage + 2 * mcse),
    linewidth = 0.5, alpha = 0.6,
    position  = position_dodge(width = 0.3)
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  scale_x_discrete(labels = overlap_x_labels) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0.7, 1.0, by = 0.1)
  ) +
  scale_colour_manual(values = estimator_colours, name = "Estimator") +
  labs(
    x     = "Degree of overlap",
    y     = "Coverage"
#    title = "Coverage of 95% CIs across overlap, error structure, and common-variance scenarios"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.text      = element_text(size = 9)
  ) +
  coord_flip( ylim = c(0.7, 1.05 ) )

print( p_coverage )


fig_path <- file.path(fig_dir, "coverage_by_overlap_het.pdf")
ggsave(fig_path,
       p_coverage,
       width = 8, height = 6)
cat("Saved figure: ", fig_path, "\n")




# And no tx variation version ----

p_coverage <- ggplot(
  filter( coverage_tbl, tx_type != "het" ),
  aes(x = overlap_label, y = coverage,
      colour = method, group = method)
) +
  facet_grid( error_type ~ k_match + common_label ) +
  geom_linerange(
    aes(ymin = coverage - 2 * mcse,
        ymax = coverage + 2 * mcse),
    linewidth = 0.5, alpha = 0.6,
    position  = position_dodge(width = 0.3)
  ) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  scale_x_discrete(labels = overlap_x_labels) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0.7, 1.0, by = 0.1)
  ) +
  scale_colour_manual(values = estimator_colours, name = "Estimator") +
  labs(
    x     = "Degree of overlap",
    y     = "Coverage"
    #    title = "Coverage of 95% CIs across overlap, error structure, and common-variance scenarios"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.text      = element_text(size = 9)
  ) +
  coord_flip( ylim = c(0.7, 1.05 ) )

print( p_coverage )


fig_path <- file.path(fig_dir, "coverage_by_overlap_homo.pdf")
ggsave(fig_path,
       p_coverage,
       width = 8, height = 6)
cat("Saved figure: ", fig_path, "\n")



resH <- filter( res, tx_type =="constant", common_label=="common" )
table( round( resH$SATT, digits = 1 ) )

# stability
resH$SATT = 3

stats <- resH %>%
  filter( method != "pop_het" & method != "pop_homo" ) %>%
  group_by( overlap_label, error_type, common_label, k_match, method ) %>%
  summarise( ESEhat = sqrt( mean( SE^2 ) ),
             trueSE = sd( att_est ),
             calibration = ESEhat / trueSE,
             R = n(),
             .groups = "drop" )
stats

stats_nested <- resH %>%
  filter(method != "pop_het", method != "pop_homo") %>%
  group_by(overlap_label, error_type, common_label, k_match, method) %>%
  nest()
stats_nested

vars <- stats_nested %>%
  mutate( vstat = map( data, ~ calc_relative_var( .x, att_est, SE^2,
                                                  criteria = "relative bias" ) ) ) %>%
  dplyr::select( -data ) %>%
  unnest( cols=vstat )
vars

stats <- left_join( stats, vars )
stats <- stats %>%
  mutate( calibration = rel_bias_var,
          mcse = rel_bias_var_mcse )


p_calibration <- ggplot(
  stats,
  aes(x = overlap_label, y = calibration,
      colour = method, group = method) ) +
  facet_grid( error_type ~ k_match + common_label ) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  geom_linerange(
    aes(ymin = calibration - 2 * mcse,
        ymax = calibration + 2 * mcse),
    linewidth = 0.5, alpha = 0.6,
    position  = position_dodge(width = 0.3)
  ) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  scale_x_discrete(labels = overlap_x_labels) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
  ) +
  scale_colour_manual(values = estimator_colours, name = "Estimator") +
  labs(
    x     = "Degree of overlap",
    y     = "Calibration"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.text      = element_text(size = 9)
  ) +
  coord_flip()

print( p_calibration )
fig_path <- file.path(fig_dir, "calibration_by_overlap.pdf")
ggsave(fig_path,
       p_calibration,
       width = 8, height = 6)
cat("Saved figure: ", fig_path, "\n")





