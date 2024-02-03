# fig_toy_example.R
# Purpose: Illustrate the importance of joint balance
#   using a toy example.
# Paper location:
#   Section 2
# Dependencies:
#   External Packages: tidyverse, latex2exp, mvtnorm
#   Internal Scripts: analysis/plot_toy.R

require(ggplot2)
# require(mvtnorm)

source("analysis/plot_toy.R")

# Load background ----------
toy_background_df_path <- "data/toy_background_df.rds"
background_df <- tryCatch({
  readRDS(toy_background_df_path)
}, error = function(err){
  message("The data file for background_df is not found. Generating one...")
  source("analysis/generate_toy_background_df.R")
  readRDS(toy_background_df_path)
}
)
background_plot <-
  create_background_plot(background_df = background_df)


## set global constants for toy examples 1, 2 --------
X1 <- 0.75; X2 <- 2.25
Y1 <- 0.75; Y2 <- 2.25
NUDGE <- 0.125

### toy example 1 -----------------------------------------------------------
## Make data
tx_dat <- tibble(
  x1 = c(X1, X2),
  x2 = c(Y1, Y2),
  z  = 1
)
co_dat <- tibble(
  x1 = c(X1, X2),
  x2 = c(Y2, Y1),
  z  = 0
)
four_points_df <- rbind(tx_dat, co_dat)

source("analysis/plot_toy.R")
toy_example_four_points_plot <-
  create_toy_example_four_points_plot(
    background_plot, four_points_df)
ggsave("writeup/figures/toyexample1.png", width=3, height=3)


# toy example 2 -----------------------------------------------------------
local_controls <- tibble(
  x1 = c(X1-1.5*NUDGE, X2-1.5*NUDGE),
  x2 = c(X1-0.5*NUDGE, X2+1.5*NUDGE),
  z  = 0
)

six_points_df <-
  rbind(four_points_df, local_controls)

toy_example_six_points_plot <-
  create_toy_example_six_points_plot(
    background_plot, six_points_df)

ggsave("writeup/figures/toyexample2.png", width=3, height=3)


# plot dataset ------------------------------------------------------------

nc <- 500
nt <- 100
f0_sd <- 1
set.seed(4)
source("R/sim_data.R")
source("R/diagnostic_plots.R")
df <- gen_df_adv(
  nc=nc,
  nt=nt,
  f0_sd = f0_sd,
  tx_effect_fun = function(X1, X2) {3*X1+3*X2},
  # f0_fun = function(x,y) {abs(x-y)})
  f0_fun = function(x,y) {
    matrix(c(x,y), ncol=2) %>%
      dmvnorm(mean = c(0.5,0.5),
              sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
  })

create_toy_df_plot(df)

ggsave("writeup/figures/sim_toy_ex.png", width=4, height=3)


# conditional expectation function ----------------------------------------

# make plot
ggplot(dat, aes(x1,x2)) +

  geom_point(aes(pch=as.factor(z)),
             size=2) +

  # label points

  annotate("text", x=X1+NUDGE, y=Y1+NUDGE, label="t[1]", parse=T) +
  annotate("text", x=X2-NUDGE, y=Y2-NUDGE, label="t[2]", parse=T) +
  annotate("text", x=X1+NUDGE, y=Y2+NUDGE, label="c[1]", parse=T) +
  annotate("text", x=X2+NUDGE, y=Y1+NUDGE, label="c[2]", parse=T) +
  annotate("text", x=X1-2*NUDGE, y=Y1-2*NUDGE, label="c[3]", parse=T) +
  annotate("text", x=X2+2*NUDGE, y=Y2+2*NUDGE, label="c[4]", parse=T) +

  # axis lines
  geom_segment(aes(x=0, y=0, xend=0, yend=3),
               arrow = arrow(length=unit(0.2, "cm"))) +
  geom_segment(aes(x=0, y=0, xend=3, yend=0),
               arrow = arrow(length=unit(0.2, "cm"))) +

  # little caliper circles
  ggforce::geom_circle(aes(x0=X1, y0=Y1, r=NUDGE*4),
                       lty = "dotted") +
  ggforce::geom_circle(aes(x0=X2, y0=Y2, r=NUDGE*4),
                       lty = "dotted") +

  # dimension/axis specification
  #  - using 0.2, 0.05 slightly shifts axis to look nicer
  coord_cartesian(xlim = c(0.1,3),
                  ylim = c(0.05,3)) +
  # scale_x_continuous(breaks = 0:5,
  #                    labels = 0:5) +
  # scale_y_continuous(breaks = 0:2,
  #                    labels = 0:2) +

  # theme it!
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank()
    # legend.position = c(2.15,.72),
    # legend.background = element_rect(color="black")
  ) +
  labs(y = TeX("$X_2$"),
       x = TeX("$X_1$"),
       pch = TeX(" $Z$"))




