# scripts/figs/fig_control_weights.R
# Visualize total control weight w_j and local residual variance s_j^2
# for one DGP draw.
# DGP: sims-variance spec (make_csm_toy_df, mid overlap)
# Matching: adaptive caliper, 5nn, scaling = (2*nbins)/(range) for numeric

devtools::load_all()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mvtnorm)
  library(tibble)
  library(patchwork)
})

# ---- DGP (from 0_sim_inference_utils.R) --------------------------------

f0_fun_xy     <- function(x, y) {
  mvtnorm::dmvnorm(cbind(x, y),
                   mean  = c(0.5, 0.5),
                   sigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)) * 20
}
tx_effect_fun_xy <- function(x, y) { 3 * x + 3 * y }

f0_fun_mat  <- function(X) { X <- as.matrix(X); f0_fun_xy(X[,1], X[,2]) }
tx_fun_mat  <- function(X) { X <- as.matrix(X); tx_effect_fun_xy(X[,1], X[,2]) }

set.seed(42)
df <- gen_df_adv_k(
  nc = 500, nt = 100, k = 2,
  f0_sd          = 0.5,
  tx_effect_fun  = tx_fun_mat,
  f0_fun         = f0_fun_mat,
  ctr_dist       = 0.5,
  prop_nc_unif   = 1/3       # "mid" overlap
)
df <- df %>% mutate(Z = as.integer(Z))

# ---- Scaling (sims-variance spec) --------------------------------------

nbins <- 6
scaling <- df %>%
  summarise(across(
    starts_with("X"),
    function(x) if (is.numeric(x)) (2 * nbins) / (max(x) - min(x)) else 1000
  ))

# ---- Matching: adaptive caliper, k = 5nn -------------------------------

mtch <- get_cal_matches(
  data       = df,
  Z ~ X1 + X2,
  rad_method = "adaptive",
  scaling    = scaling,
  k          = 5,
  warn       = FALSE,
  est_method = "scm"
)

# ---- Extract total control weight w_j per control unit -----------------

co_weights <- result_table(mtch, return = "agg_co_units") %>%
  filter(Z == 0)                     # keep only control rows

# Join back covariates (id is character in the result table)
df_ctrl <- df %>%
  filter(Z == 0) %>%
  mutate(id = as.character(id))

plot_df <- df_ctrl %>%
  left_join(co_weights %>% select(id, w_j = weights), by = "id") %>%
  mutate(w_j = replace(w_j, is.na(w_j), 0),   # unmatched controls get 0
         matched = w_j > 1e-9)

treated_df <- df %>% filter(Z == 1)

# ---- Compute s_j^2 and Cov_p(w_j, s_j^2) ------------------------------

matches_full <- full_unit_table(mtch)   # subclass-level data needed for s_j^2

sj_result <- calculate_s_j_sq(matches_full, outcome = "Y", treatment = "Z")

# Join s_j^2 to plot_df
plot_df <- plot_df %>%
  left_join(sj_result$s_j_sq %>% mutate(id = as.character(id)), by = "id")

# Empirical covariance Cov_p(w_j, s_j^2) over matched controls
cov_w_s <- with(
  plot_df %>% filter(matched, !is.na(s_j_sq)),
  cov(w_j, s_j_sq)
)
message(sprintf("Cov_p(w_j, s_j^2) = %.4f", cov_w_s))

# ---- Plot 1: total control weight w_j ----------------------------------

p1 <- ggplot() +
  geom_point(
    data  = plot_df %>% filter(!matched),
    aes(X1, X2),
    color = "grey80", size = 1.5, alpha = 0.6, shape = 16
  ) +
  geom_point(
    data  = plot_df %>% filter(matched),
    aes(X1, X2, color = w_j, size = w_j),
    alpha = 0.85, shape = 16
  ) +
  geom_point(
    data  = treated_df,
    aes(X1, X2),
    color = "firebrick", size = 2.5, shape = 17, alpha = 0.9
  ) +
  scale_color_viridis_c(
    option    = "plasma",
    name      = expression(w[j]),
    direction = -1,
    trans     = "sqrt"
  ) +
  scale_size_continuous(
    name   = expression(w[j]),
    range  = c(1, 5),
    trans  = "sqrt",
    guide  = "none"
  ) +
  labs(
    title    = expression("Total control weight " * w[j]),
    subtitle = "Adaptive caliper (5nn) · mid overlap · grey = unmatched · red triangles = treated",
    x = expression(X[1]), y = expression(X[2])
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 9, color = "grey40")
  )

# ---- Plot 2: local residual variance s_j^2 -----------------------------

p2 <- ggplot() +
  geom_point(
    data  = plot_df %>% filter(!matched),
    aes(X1, X2),
    color = "grey80", size = 1.5, alpha = 0.6, shape = 16
  ) +
  geom_point(
    data  = plot_df %>% filter(matched, !is.na(s_j_sq)),
    aes(X1, X2, color = s_j_sq, size = s_j_sq),
    alpha = 0.85, shape = 16
  ) +
  geom_point(
    data  = treated_df,
    aes(X1, X2),
    color = "firebrick", size = 2.5, shape = 17, alpha = 0.9
  ) +
  scale_color_viridis_c(
    option    = "viridis",
    name      = expression(s[j]^2),
    direction = -1,
    trans     = "sqrt"
  ) +
  scale_size_continuous(
    name   = expression(s[j]^2),
    range  = c(1, 5),
    trans  = "sqrt",
    guide  = "none"
  ) +
  labs(
    title    = expression("Local residual variance " * s[j]^2),
    subtitle = sprintf("Cov_p(w[j], s[j]^2) = %.4f", cov_w_s),
    x = expression(X[1]), y = expression(X[2])
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 9, color = "grey40")
  )

# ---- Combine and save (homoskedastic baseline) -------------------------

combined <- p1 + p2 + plot_layout(ncol = 2)
print(combined)

out_path <- here::here("scripts/figs/fig_control_weights.pdf")
ggsave(out_path, combined, width = 14, height = 5.5)
message("Saved: ", out_path)

# ========================================================================
# Het-sigma scenarios: σ₀(x)=0.2+(x1-x2)²; scenario 2 adds σ₁=σ₀+0.5
# ========================================================================

source(here::here("scripts/sims-variance-het-sigma/0_sim_inference_utils.R"))

make_het_panels <- function(sigma1_extra, seed_val = 42, scenario_label = "") {
  set.seed(seed_val)
  df_h <- make_het_sigma_df(
    nc = 500, nt = 100,
    prop_nc_unif = 1/3,   # mid overlap
    ctr_dist     = 0.5,
    seed         = seed_val,
    sigma1_extra = sigma1_extra
  )

  scaling_h <- compute_toy_scaling(df_h)

  mtch_h <- get_cal_matches(
    data       = df_h,
    Z ~ X1 + X2,
    rad_method = "adaptive",
    scaling    = scaling_h,
    k          = 2,
    warn       = FALSE,
    est_method = "scm"
  )

  co_weights_h <- result_table(mtch_h, return = "agg_co_units") %>%
    filter(Z == 0)

  df_ctrl_h <- df_h %>%
    filter(Z == 0) %>%
    mutate(id = as.character(id))

  plot_df_h <- df_ctrl_h %>%
    left_join(co_weights_h %>% select(id, w_j = weights), by = "id") %>%
    mutate(
      w_j     = replace(w_j, is.na(w_j), 0),
      matched = w_j > 1e-9
    )

  treated_h <- df_h %>% filter(Z == 1)

  matches_full_h <- full_unit_table(mtch_h)
  sj_h <- calculate_s_j_sq(matches_full_h, outcome = "Y", treatment = "Z")

  plot_df_h <- plot_df_h %>%
    left_join(sj_h$s_j_sq %>% mutate(id = as.character(id)), by = "id")

  cov_ws_h <- with(
    plot_df_h %>% filter(matched, !is.na(s_j_sq)),
    cov(w_j, s_j_sq)
  )
  message(sprintf("[%s] Cov_p(w_j, s_j^2) = %.4f", scenario_label, cov_ws_h))

  subtitle_base <- sprintf("%s · adaptive (2nn) · mid overlap · grey = unmatched", scenario_label)

  ph1 <- ggplot() +
    geom_point(data = plot_df_h %>% filter(!matched),
               aes(X1, X2), color = "grey80", size = 1.5, alpha = 0.6, shape = 16) +
    geom_point(data = plot_df_h %>% filter(matched),
               aes(X1, X2, color = w_j, size = w_j),
               alpha = 0.85, shape = 16) +
    geom_point(data = treated_h,
               aes(X1, X2), color = "firebrick", size = 2.5, shape = 17, alpha = 0.9) +
    scale_color_viridis_c(option = "plasma", name = expression(w[j]),
                          direction = -1, trans = "sqrt") +
    scale_size_continuous(name = expression(w[j]), range = c(1, 5),
                          trans = "sqrt", guide = "none") +
    labs(title    = expression("Control weight " * w[j]),
         subtitle = subtitle_base,
         x = expression(X[1]), y = expression(X[2])) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 8, color = "grey40"))

  ph2 <- ggplot() +
    geom_point(data = plot_df_h %>% filter(!matched),
               aes(X1, X2), color = "grey80", size = 1.5, alpha = 0.6, shape = 16) +
    geom_point(data = plot_df_h %>% filter(matched, !is.na(s_j_sq)),
               aes(X1, X2, color = s_j_sq, size = s_j_sq),
               alpha = 0.85, shape = 16) +
    geom_point(data = treated_h,
               aes(X1, X2), color = "firebrick", size = 2.5, shape = 17, alpha = 0.9) +
    scale_color_viridis_c(option = "viridis", name = expression(s[j]^2),
                          direction = -1, trans = "sqrt") +
    scale_size_continuous(name = expression(s[j]^2), range = c(1, 5),
                          trans = "sqrt", guide = "none") +
    labs(title    = expression("Local variance " * s[j]^2),
         subtitle = sprintf("Cov_p(w[j], s[j]^2) = %.4f", cov_ws_h),
         x = expression(X[1]), y = expression(X[2])) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 8, color = "grey40"))

  list(p1 = ph1, p2 = ph2)
}

panels_s1 <- make_het_panels(
  sigma1_extra   = 0,
  seed_val       = 42,
  scenario_label = "Scenario 1: sigma0=sigma1=0.2+(x1-x2)^2"
)

panels_s2 <- make_het_panels(
  sigma1_extra   = 0.5,
  seed_val       = 42,
  scenario_label = "Scenario 2: sigma0=0.2+(x1-x2)^2, sigma1=sigma0+0.5"
)

het_combined <- (panels_s1$p1 + panels_s1$p2) /
                (panels_s2$p1 + panels_s2$p2) +
  plot_layout(ncol = 1)

print(het_combined)

out_het <- here::here("scripts/figs/fig_het_sigma_weights.pdf")
ggsave(out_het, het_combined, width = 14, height = 11)
message("Saved: ", out_het)
