#!/usr/bin/env Rscript
# Illustration of "fair" vs "current" parameterization when comparing CSM to CEM.
#
# The proposition in §sec:compCEM states:
#   CSM with caliper (radius) = π/2  and  CEM with bin-size = π
#   have the same catchment area, but CSM has half the bias bound.
#
# In both panels we set:
#   CEM  bin-width = R / nbins
#   CSM  caliper (radius) = R / nbins      ← same numeric value
#
# Panel (a) — CURRENT CODE (same nbins for both):
#   CEM catchment width = R/nbins (one bin)
#   CSM catchment width = 2·R/nbins (radius on each side)  → CSM is 2× larger
#
# Panel (b) — FAIR COMPARISON (CEM uses nbins/2):
#   CEM bin-width = 2·R/nbins  → catchment width = 2·R/nbins
#   CSM caliper   =   R/nbins  → catchment width = 2·R/nbins  (equal!)
#   But CSM max dist = R/nbins,  CEM max dist = 2·R/nbins  → CSM bias bound 2× tighter

library(tidyverse)
library(patchwork)
library(latex2exp)

set.seed(42)

# ── Setup ─────────────────────────────────────────────────────────────────────
RANGE <- 2        # covariate range [0, RANGE]
NBINS <- 4        # nbins used for CSM in the current simulation code
r_csm <- RANGE / NBINS   # CSM caliper (radius) = R/nbins = 0.5

# Treated unit — placed near a bin edge so the CEM boundary effect is visible
tx <- c(x1 = 0.75, x2 = 0.75)

# ── Random control points ──────────────────────────────────────────────────────
n_co <- 60
pts <- tibble(
  x1 = runif(n_co, 0, RANGE),
  x2 = runif(n_co, 0, RANGE),
  z  = 0L
) %>%
  bind_rows(tibble(x1 = tx["x1"], x2 = tx["x2"], z = 1L))

# ── Helpers ────────────────────────────────────────────────────────────────────
cem_bin_lo <- function(v, width) floor(v / width) * width

cem_matched <- function(px, py, tx, ty, width) {
  cem_bin_lo(px, width) == cem_bin_lo(tx["x1"], width) &
    cem_bin_lo(py, width) == cem_bin_lo(tx["x2"], width)
}
csm_matched <- function(px, py, radius) {
  abs(px - tx["x1"]) <= radius & abs(py - tx["x2"]) <= radius
}

# ── Panel builder ──────────────────────────────────────────────────────────────
# cem_nbins  : number of CEM bins (controls bin-width = RANGE/cem_nbins)
# csm_radius : CSM caliper — always R/NBINS = 0.5
make_panel <- function(cem_nbins, csm_radius, panel_label,
                       cem_label, csm_label, catchment_note) {

  bw <- RANGE / cem_nbins   # CEM bin-width

  pts_aug <- pts %>%
    filter(z == 0) %>%
    mutate(
      in_cem = cem_matched(x1, x2, tx, tx, bw),
      in_csm = csm_matched(x1, x2, csm_radius),
      match_status = case_when(
        in_cem & in_csm ~ "Both",
        in_cem          ~ "CEM only",
        in_csm          ~ "CSM only",
        TRUE            ~ "Neither"
      )
    )

  grid_x <- seq(0, RANGE, by = bw)
  grid_y <- seq(0, RANGE, by = bw)

  # CEM bin containing treated unit
  bx_lo <- cem_bin_lo(tx["x1"], bw)
  by_lo <- cem_bin_lo(tx["x2"], bw)

  # CSM caliper box
  cx_lo <- tx["x1"] - csm_radius;  cx_hi <- tx["x1"] + csm_radius
  cy_lo <- tx["x2"] - csm_radius;  cy_hi <- tx["x2"] + csm_radius

  # CEM max-distance arrow: from treated unit to far edge of its bin
  cem_arrow_y <- by_lo + bw * 0.25   # position arrow in lower quarter of bin

  # CSM radius arrow: from treated unit rightward to caliper edge
  csm_arrow_y <- tx["x2"] - csm_radius * 0.5

  status_colors <- c(
    "Both"     = "#8B008B",
    "CEM only" = "#E69F00",
    "CSM only" = "#009E73",
    "Neither"  = "grey65"
  )
  status_shapes <- c(
    "Both"     = 16,
    "CEM only" = 16,
    "CSM only" = 16,
    "Neither"  = 1
  )

  ggplot() +
    # ── Catchment shading ──────────────────────────────────────────────────
    annotate("rect",
             xmin = bx_lo, xmax = bx_lo + bw,
             ymin = by_lo, ymax = by_lo + bw,
             fill = "#E69F00", alpha = 0.18, color = NA) +
    annotate("rect",
             xmin = cx_lo, xmax = cx_hi,
             ymin = cy_lo, ymax = cy_hi,
             fill = "#009E73", alpha = 0.18, color = NA) +
    # ── CEM grid ──────────────────────────────────────────────────────────
    geom_vline(xintercept = grid_x, linetype = "dashed",
               color = "#C07000", linewidth = 0.55) +
    geom_hline(yintercept = grid_y, linetype = "dashed",
               color = "#C07000", linewidth = 0.55) +
    # ── CSM caliper border ────────────────────────────────────────────────
    annotate("rect",
             xmin = cx_lo, xmax = cx_hi,
             ymin = cy_lo, ymax = cy_hi,
             fill = NA, color = "#009E73", linewidth = 1.3) +
    # ── Points ────────────────────────────────────────────────────────────
    geom_point(data = pts_aug,
               aes(x1, x2, color = match_status, shape = match_status),
               size = 2.4) +
    geom_point(data = tibble(x1 = tx["x1"], x2 = tx["x2"]),
               aes(x1, x2), shape = 17, color = "red", size = 4.5) +
    annotate("text", x = tx["x1"] + 0.07, y = tx["x2"] + 0.10,
             label = expression(italic(t)[1]), size = 3.8, fontface = "italic") +
    # ── CEM max-distance arrow: treated → far edge of bin ─────────────────
    annotate("segment",
             x = tx["x1"], xend = bx_lo + bw,
             y = cem_arrow_y, yend = cem_arrow_y,
             arrow = arrow(length = unit(0.18, "cm"), ends = "both"),
             color = "#C07000", linewidth = 1.0) +
    annotate("text",
             x = (tx["x1"] + bx_lo + bw) / 2,
             y = cem_arrow_y + 0.09,
             label = cem_label,
             color = "#C07000", size = 3.0, fontface = "bold") +
    # ── CSM caliper arrow: treated → right edge ────────────────────────────
    annotate("segment",
             x = tx["x1"], xend = cx_hi,
             y = csm_arrow_y, yend = csm_arrow_y,
             arrow = arrow(length = unit(0.18, "cm"), ends = "last"),
             color = "#007050", linewidth = 1.0) +
    annotate("text",
             x = (tx["x1"] + cx_hi) / 2,
             y = csm_arrow_y - 0.10,
             label = csm_label,
             color = "#007050", size = 3.0, fontface = "bold") +
    # ── Catchment-comparison note ──────────────────────────────────────────
    annotate("text",
             x = RANGE / 2, y = RANGE - 0.08,
             label = catchment_note,
             size = 3.1, hjust = 0.5, vjust = 1, color = "grey20") +
    # ── Scales / theme ────────────────────────────────────────────────────
    scale_color_manual(values = status_colors, name = "Matched by") +
    scale_shape_manual(values = status_shapes, name = "Matched by") +
    coord_fixed(xlim = c(0, RANGE), ylim = c(0, RANGE)) +
    scale_x_continuous(breaks = seq(0, RANGE, 0.5),
                       labels = function(x) ifelse(x %% 1 == 0, x, "")) +
    scale_y_continuous(breaks = seq(0, RANGE, 0.5),
                       labels = function(x) ifelse(x %% 1 == 0, x, "")) +
    labs(
      title = panel_label,
      x = TeX("$X_1$"),
      y = TeX("$X_2$")
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "bottom",
      legend.text     = element_text(size = 8),
      legend.title    = element_text(size = 9),
      plot.title      = element_text(face = "bold", size = 11)
    )
}

# ── Panel A: current code — same nbins for CEM and CSM ────────────────────────
# Both: CEM bin-width = R/nbins = 0.5,  CSM caliper = R/nbins = 0.5
pA <- make_panel(
  cem_nbins      = NBINS,
  csm_radius     = r_csm,
  panel_label    = sprintf("(a)  Current code:  same nbins = %d for both", NBINS),
  cem_label      = "CEM max dist  =  R/nbins",
  csm_label      = "CSM caliper  =  R/nbins",
  catchment_note = paste0(
    "CEM bin-width = R/nbins = ", r_csm, "  →  CEM catchment width = R/nbins\n",
    "CSM caliper (radius) = R/nbins = ", r_csm, "  →  CSM catchment width = 2×R/nbins\n",
    "CSM catchment is 2× larger  |  bias bounds are EQUAL"
  )
)

# ── Panel B: fair comparison — CEM uses nbins/2, CSM stays at nbins ───────────
# CEM bin-width = 2·R/nbins = 1.0,  CSM caliper = R/nbins = 0.5
pB <- make_panel(
  cem_nbins      = NBINS / 2,
  csm_radius     = r_csm,
  panel_label    = sprintf("(b)  Fair comparison:  CEM uses nbins/2 = %d,  CSM keeps nbins = %d",
                           NBINS / 2, NBINS),
  cem_label      = "CEM max dist  =  2×R/nbins",
  csm_label      = "CSM caliper  =  R/nbins",
  catchment_note = paste0(
    "CEM bin-width = 2×R/nbins = ", 2 * r_csm, "  →  CEM catchment width = 2×R/nbins\n",
    "CSM caliper (radius) = R/nbins = ", r_csm, "  →  CSM catchment width = 2×R/nbins\n",
    "Catchments are EQUAL  |  CSM bias bound is 2× tighter (Proposition §sec:compCEM)"
  )
)

fig <- pA + pB +
  plot_annotation(
    title   = "CSM vs CEM: catchment areas under equal-nbins vs fair parameterization",
    caption = paste(
      "Red triangle: treated unit t1.  Orange dashed: CEM bin boundaries; orange shading: CEM catchment.",
      "Green box: CSM L∞ caliper; green shading: CSM catchment.  nbins = 4, R = 2.",
      sep = "  "
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 12),
      plot.caption = element_text(size = 7.5, hjust = 0)
    )
  )

out_dir <- here::here("writeup/figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "cem_csm_fair_comparison.pdf"),
       fig, width = 13, height = 6.5)
ggsave(file.path(out_dir, "cem_csm_fair_comparison.png"),
       fig, width = 13, height = 6.5, dpi = 150)
message("Saved to ", out_dir)
