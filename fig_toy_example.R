
# illustrate toy example of why joint balance is important

require(tidyverse)
require(latex2exp)


# set global constants

X1 <- 0.75
X2 <- 2.25
Y1 <- 0.75
Y2 <- 2.25
NUDGE <- 0.125



# build gradient-looking background ---------------------------------------

require(mvtnorm)
MU <- c(1.5,1.5)
SIG <- matrix(c(1,0.5,0.5,1), nrow=2)

dat_background <- expand.grid(
  x1 = seq(0.025,3,by=0.01),
  x2 = seq(0.025,3,by=0.01)
  # z  = abs(x1-x2)    # straight diagonal gradients
) %>% 
  rowwise() %>% 
  mutate(z = dmvnorm(x = c(x1,x2),
                     mean = MU,
                     sigma = SIG))


# toy example 1 -----------------------------------------------------------

# make fake data

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
dat1 <- rbind(tx_dat, co_dat)


# make plot
ggplot(dat1, aes(x1,x2)) +
  
  # contours
  # TODO: maybe roundish contours are clearer?
  geom_tile(data=dat_background,
            aes(fill=z),
            alpha=0.65) +
  scale_fill_gradient(low="blue", high="orange") +
  geom_contour(data=dat_background,
               aes(z=z),
               color="black",
               alpha = 0.1,
               lty = "solid",
               bins = 5) +
  
  geom_point(aes(pch=as.factor(z)),
             size=2) +
  
  # label points
  annotate("text", x=X1+NUDGE, y=Y1+NUDGE, label="t[1]", parse=T) +
  annotate("text", x=X2+NUDGE, y=Y2+NUDGE, label="t[2]", parse=T) +
  annotate("text", x=X1+NUDGE, y=Y2+NUDGE, label="c[1]", parse=T) +
  annotate("text", x=X2+NUDGE, y=Y1+NUDGE, label="c[2]", parse=T) +
  
  # axis lines
  geom_segment(aes(x=0, y=0, xend=0, yend=3),
               arrow = arrow(length=unit(0.2, "cm")),
               linewidth=1) +
  geom_segment(aes(x=0, y=0, xend=3, yend=0),
               arrow = arrow(length=unit(0.2, "cm")),
               linewidth=1) +
  
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

ggsave("writeup/figures/toyexample1.png", width=3, height=3)




# toy example 2 -----------------------------------------------------------

# make fake data
tx_dat <- tibble(
  x1 = c(X1, X2),
  x2 = c(Y1, Y2),
  z  = 1
)
co_dat <- tibble(
  x1 = c(X1, X2, X1-1.5*NUDGE, X2-1.5*NUDGE),
  x2 = c(Y2, Y1, X1-0.5*NUDGE, X2+1.5*NUDGE),
  z  = 0
)
dat2 <- rbind(tx_dat, co_dat)


# make plot
ggplot(dat2, aes(x1,x2)) +

  # contours
  # TODO: maybe roundish contours are clearer?
  geom_tile(data=dat_background,
            aes(fill=z),
            alpha=0.65) +
  scale_fill_gradient(low="blue", high="orange") +
  geom_contour(data=dat_background,
               aes(z=z),
               color="black",
               alpha = 0.1,
               lty = "solid",
               bins = 5) +
  
  geom_point(aes(pch=as.factor(z)),
             size=2) +
  
  # label points
  
  annotate("text", x=X1+NUDGE, y=Y1+NUDGE, label="t[1]", parse=T) +
  annotate("text", x=X2+NUDGE, y=Y2+NUDGE, label="t[2]", parse=T) +
  annotate("text", x=X1+NUDGE, y=Y2+NUDGE, label="c[1]", parse=T) +
  annotate("text", x=X2+NUDGE, y=Y1+NUDGE, label="c[2]", parse=T) +
  annotate("text", x=X1-2.9*NUDGE, y=Y1-0.5*NUDGE, label="c[3]", parse=T) +
  annotate("text", x=X2-2.3*NUDGE, y=Y2+0.7*NUDGE, label="c[4]", parse=T) +
  
  # axis lines
  geom_segment(aes(x=0, y=0, xend=0, yend=3),
               arrow = arrow(length=unit(0.2, "cm")),
               linewidth=1) +
  geom_segment(aes(x=0, y=0, xend=3, yend=0),
               arrow = arrow(length=unit(0.2, "cm")),
               linewidth=1) +
  
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

ggsave("writeup/figures/toyexample2.png", width=3, height=3)



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




