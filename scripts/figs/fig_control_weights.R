# scripts/figs/fig_control_weights.R
# Visualize total control weight w_j for one DGP draw
# DGP: sims-variance spec (make_csm_toy_df, mid overlap)
# Matching: adaptive caliper, 5nn, scaling = (2*nbins)/(range) for numeric

devtools::load_all()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mvtnorm)
  library(tibble)
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

# ---- Plot --------------------------------------------------------------

p <- ggplot() +
  # unmatched controls (w_j == 0) in light grey
  geom_point(
    data  = plot_df %>% filter(!matched),
    aes(X1, X2),
    color = "grey80", size = 1.5, alpha = 0.6, shape = 16
  ) +
  # matched controls, coloured by w_j
  geom_point(
    data  = plot_df %>% filter(matched),
    aes(X1, X2, color = w_j, size = w_j),
    alpha = 0.85, shape = 16
  ) +
  # treated units
  geom_point(
    data  = treated_df,
    aes(X1, X2),
    color = "firebrick", size = 2.5, shape = 17, alpha = 0.9
  ) +
  scale_color_viridis_c(
    option    = "plasma",
    name      = expression(w[j]),
    direction = -1,
    trans     = "sqrt"         # sqrt-scale helps reveal low-weight controls
  ) +
  scale_size_continuous(
    name   = expression(w[j]),
    range  = c(1, 5),
    trans  = "sqrt",
    guide  = "none"
  ) +
  labs(
    title    = expression("Total control weight " * w[j]),
    subtitle = "Adaptive caliper (5nn) · mid overlap · grey = unmatched controls · red triangles = treated",
    x = expression(X[1]), y = expression(X[2])
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 9, color = "grey40")
  )

print(p)

out_path <- here::here("scripts/figs/fig_control_weights.pdf")
ggsave(out_path, p, width = 7, height = 5.5)
message("Saved: ", out_path)
