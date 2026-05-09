# scripts/figs/fig_het_sigma_bimodal.R
#
# Compares two heteroskedastic σ₀ designs:
#   Design A (diagonal):  σ₀(x) = 0.2 + (x1-x2)²
#   Design B (bimodal):   σ₀(x) = 0.2 + exp(-‖x-(0.25,0.25)‖²/0.08)
#                                      + exp(-‖x-(0.75,0.75)‖²/0.08)
#
# For each design: σ₀ contour | w_j scatter | s_j² scatter
# (6 panels total, 2 rows × 3 cols)

devtools::load_all()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(patchwork)
  library(here)
  library(mvtnorm)
})

source(here::here("scripts/sims-variance-het-sigma/0_sim_inference_utils.R"))
source(here::here("scripts/sims-variance-het-sigma-bimodal/0_sim_inference_utils.R"))

# -------------------------------------------------------------------
# Helper: make a σ₀ contour plot for either design
# -------------------------------------------------------------------

make_sigma_contour <- function(sigma_fun, title_str,
                               xlim = c(0, 1), ylim = c(0, 1),
                               n_grid = 150) {
  grid_df <- expand.grid(
    X1 = seq(xlim[1], xlim[2], length.out = n_grid),
    X2 = seq(ylim[1], ylim[2], length.out = n_grid)
  ) %>%
    mutate(sigma0 = sigma_fun(X1, X2))

  ggplot(grid_df, aes(X1, X2)) +
    geom_raster(aes(fill = sigma0), interpolate = TRUE) +
    geom_contour(aes(z = sigma0), color = "white", alpha = 0.55,
                 linewidth = 0.35) +
    # mark the two treated cluster centres
    annotate("point", x = c(0.25, 0.75), y = c(0.25, 0.75),
             shape = 17, size = 3, color = "firebrick") +
    annotate("point", x = c(0.75, 0.25), y = c(0.25, 0.75),
             shape = 16, size = 2.5, color = "white", alpha = 0.7) +
    scale_fill_viridis_c(name = expression(sigma[0](x)), option = "magma") +
    labs(title = title_str,
         subtitle = "▲ treated centres  ● control centres",
         x = expression(X[1]), y = expression(X[2])) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 8, color = "grey40"))
}

# -------------------------------------------------------------------
# Helper: make w_j + s_j² panels from matching results
# -------------------------------------------------------------------

make_wj_sj_panels <- function(df, scaling, sigma1_extra = 0,
                               design_label = "") {
  mtch <- get_cal_matches(
    data       = df,
    Z ~ X1 + X2,
    rad_method = "adaptive",
    scaling    = scaling,
    k          = 2,
    warn       = FALSE,
    est_method = "scm"
  )

  co_wts <- result_table(mtch, return = "agg_co_units") %>% filter(Z == 0)
  df_ctrl <- df %>% filter(Z == 0) %>% mutate(id = as.character(id))

  plot_df <- df_ctrl %>%
    left_join(co_wts %>% select(id, w_j = weights), by = "id") %>%
    mutate(w_j = replace(w_j, is.na(w_j), 0), matched = w_j > 1e-9)

  treated_pts <- df %>% filter(Z == 1)

  mf <- full_unit_table(mtch)
  sj <- calculate_s_j_sq(mf, outcome = "Y", treatment = "Z")
  plot_df <- plot_df %>%
    left_join(sj$s_j_sq %>% mutate(id = as.character(id)), by = "id")

  cov_ws <- with(plot_df %>% filter(matched, !is.na(s_j_sq)),
                 cov(w_j, s_j_sq))
  message(sprintf("[%s, sigma1_extra=%.1f] Cov_p(w_j, s_j^2) = %.4f",
                  design_label, sigma1_extra, cov_ws))

  subtitle_w  <- sprintf("%s · adaptive(2nn) · mid overlap", design_label)
  subtitle_sj <- sprintf("Cov_p(w[j], s[j]^2) = %.4f", cov_ws)

  p_wj <- ggplot() +
    geom_point(data = plot_df %>% filter(!matched),
               aes(X1, X2), color = "grey80", size = 1.2, alpha = 0.5, shape = 16) +
    geom_point(data = plot_df %>% filter(matched),
               aes(X1, X2, color = w_j, size = w_j), alpha = 0.85, shape = 16) +
    geom_point(data = treated_pts,
               aes(X1, X2), color = "firebrick", size = 2.2, shape = 17, alpha = 0.9) +
    scale_color_viridis_c(option = "plasma", name = expression(w[j]),
                          direction = -1, trans = "sqrt") +
    scale_size_continuous(name = expression(w[j]), range = c(1, 4.5),
                          trans = "sqrt", guide = "none") +
    labs(title = expression("Control weight " * w[j]),
         subtitle = subtitle_w, x = expression(X[1]), y = expression(X[2])) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 8, color = "grey40"))

  p_sj <- ggplot() +
    geom_point(data = plot_df %>% filter(!matched),
               aes(X1, X2), color = "grey80", size = 1.2, alpha = 0.5, shape = 16) +
    geom_point(data = plot_df %>% filter(matched, !is.na(s_j_sq)),
               aes(X1, X2, color = s_j_sq, size = s_j_sq), alpha = 0.85, shape = 16) +
    geom_point(data = treated_pts,
               aes(X1, X2), color = "firebrick", size = 2.2, shape = 17, alpha = 0.9) +
    scale_color_viridis_c(option = "viridis", name = expression(s[j]^2),
                          direction = -1, trans = "sqrt") +
    scale_size_continuous(name = expression(s[j]^2), range = c(1, 4.5),
                          trans = "sqrt", guide = "none") +
    labs(title = expression("Local variance " * s[j]^2),
         subtitle = subtitle_sj, x = expression(X[1]), y = expression(X[2])) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 8, color = "grey40"))

  list(p_wj = p_wj, p_sj = p_sj, cov_ws = cov_ws)
}

# -------------------------------------------------------------------
# Generate data for both designs (scenario 1 only for visualization)
# -------------------------------------------------------------------

set.seed(42)
df_diag <- make_het_sigma_df(
  nc = 500, nt = 100, prop_nc_unif = 1/3, ctr_dist = 0.5,
  seed = 42, sigma1_extra = 0
)

set.seed(42)
df_bim <- make_bimodal_sigma_df(
  nc = 500, nt = 100, prop_nc_unif = 1/3, ctr_dist = 0.5,
  seed = 42, sigma1_extra = 0
)

scaling_diag <- compute_toy_scaling(df_diag)
scaling_bim  <- compute_toy_scaling(df_bim)

# -------------------------------------------------------------------
# σ₀ contour plots
# -------------------------------------------------------------------

# wrap into positional functions for make_sigma_contour
sigma_diag_fn <- function(x1, x2) { 0.2 + (x1 - x2)^2 }
sigma_bim_fn  <- function(x1, x2) { sigma0_bimodal(x1, x2) }

p_sigma_diag <- make_sigma_contour(
  sigma_diag_fn,
  title_str = expression("Design A: " * sigma[0](x) == 0.2 + (x[1] - x[2])^2)
)

p_sigma_bim <- make_sigma_contour(
  sigma_bim_fn,
  title_str = expression("Design B: " * sigma[0](x) == "bimodal Gaussian bumps")
)

# -------------------------------------------------------------------
# w_j and s_j² panels
# -------------------------------------------------------------------

panels_diag <- make_wj_sj_panels(df_diag, scaling_diag,
                                  design_label = "Design A (diagonal)")
panels_bim  <- make_wj_sj_panels(df_bim,  scaling_bim,
                                  design_label = "Design B (bimodal)")

# -------------------------------------------------------------------
# Compose 2×3 figure
# -------------------------------------------------------------------

row_A <- p_sigma_diag | panels_diag$p_wj | panels_diag$p_sj
row_B <- p_sigma_bim  | panels_bim$p_wj  | panels_bim$p_sj

fig_combined <- (row_A / row_B) +
  plot_annotation(
    title   = expression("Effect of σ₀(x) design on Cov"[p] * "(w"[j] * ", s"[j]^2 * ")"),
    caption = sprintf(
      "Design A  Cov_p = %.4f   |   Design B  Cov_p = %.4f",
      panels_diag$cov_ws, panels_bim$cov_ws
    ),
    theme = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 10, color = "grey30")
    )
  )

print(fig_combined)

out_path <- here::here("scripts/figs/fig_het_sigma_bimodal.pdf")
ggsave(out_path, fig_combined, width = 18, height = 11)
message("Saved: ", out_path)
