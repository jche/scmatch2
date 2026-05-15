
# create matching example figure


require(tidyverse)
require(latex2exp)

set.seed(90210)


# make fake data
tx_dat <- tibble(
  x1 = 6-c(1.1, 1.8, 2.1, 2.2),
  x2 = c(1.5, 1.1, 0.8, 2.3),
  z  = 1
)
co_dat <- tibble(
  x1 = 6-c(1.3,2.1,2.9,3.1,3.2,3.3,3.3,3.8,4.2,4.2,4.4,4.5,4.7,4.9),
  x2 = c(1.5,1.1,0.6,2.3,0.3,1.2,0.6,2.1,1.5,2.3,1.2,0.7,2.3,0.8),
  z  = 0
)
dat <- rbind(tx_dat, co_dat)


# make plot
ggplot(dat, aes(x1,x2)) +

  # caliper
  geom_rect(xmin=2.2-1, xmax=2.2+1, ymin=1.9-0.5, ymax=1.9+0.5,
            fill = "red",
            color = "black",
            linetype = "dotted",
            alpha = 0.01) +

  geom_point(aes(pch=as.factor(z)),
             size=2) +

  # special point
  geom_point(data=tibble(x1=2.2, x2=1.9),
             pch=17, color="red", size=2) +
  annotate("text", x=2.35, y=1.9, label="t[1]", parse=T) +

  # label some control points
  annotate("text", x=2.35, y=2.1, label="c[1]", parse=T) +
  annotate("text", x=3.75, y=1.1, label="c[2]", parse=T) +

  # axis lines
  geom_segment(aes(x=0, y=0, xend=0, yend=2.5),
               arrow = arrow(length=unit(0.2, "cm"))) +
  geom_segment(aes(x=0, y=0, xend=5, yend=0),
               arrow = arrow(length=unit(0.2, "cm"))) +

  # gridlines
  geom_segment(aes(x=2, y=0, xend=2, yend=2.5),
               lty = "dashed") +
  geom_segment(aes(x=4, y=0, xend=4, yend=2.5),
               lty = "dashed") +
  geom_segment(aes(x=0, y=1, xend=5, yend=1),
               lty = "dashed") +
  geom_segment(aes(x=0, y=2, xend=5, yend=2),
               lty = "dashed") +

  # dimension/axis specification
  #  - using 0.2, 0.05 slightly shifts axis to look nicer
  coord_cartesian(xlim = c(0.2,5),
                  ylim = c(0.05,2.5)) +
  scale_x_continuous(breaks = 0:5,
                     labels = 0:5) +
  scale_y_continuous(breaks = 0:2,
                     labels = 0:2) +

  # theme it!
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = c(.1,.72),
    legend.background = element_rect(color="black")) +
  labs(y = TeX("$X_2$"),
       x = TeX("$X_1$"),
       pch = TeX(" $Z$"))

ggsave("figures/show_cem_calipers.png",
       width=6,
       height=6)





# without caliper
ggplot(dat, aes(x1,x2)) +

  geom_point(aes(pch=as.factor(z)),
             size=2) +

  # special point
  geom_point(data=tibble(x1=2.2, x2=1.9),
             pch=17, color="red", size=2) +
  annotate("text", x=2.35, y=1.9, label="t[1]", parse=T) +

  # label some control points
  annotate("text", x=2.35, y=2.1, label="c[1]", parse=T) +
  annotate("text", x=3.75, y=1.1, label="c[2]", parse=T) +

  # axis lines
  geom_segment(aes(x=0, y=0, xend=0, yend=2.5),
               arrow = arrow(length=unit(0.2, "cm"))) +
  geom_segment(aes(x=0, y=0, xend=5, yend=0),
               arrow = arrow(length=unit(0.2, "cm"))) +

  # gridlines
  geom_segment(aes(x=2, y=0, xend=2, yend=2.5),
               lty = "dashed") +
  geom_segment(aes(x=4, y=0, xend=4, yend=2.5),
               lty = "dashed") +
  geom_segment(aes(x=0, y=1, xend=5, yend=1),
               lty = "dashed") +
  geom_segment(aes(x=0, y=2, xend=5, yend=2),
               lty = "dashed") +

  # dimension/axis specification
  #  - using 0.2, 0.05 slightly shifts axis to look nicer
  coord_cartesian(xlim = c(0.2,5),
                  ylim = c(0.05,2.5)) +
  scale_x_continuous(breaks = 0:5,
                     labels = 0:5) +
  scale_y_continuous(breaks = 0:2,
                     labels = 0:2) +

  # theme it!
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = c(.1,.72),
    legend.background = element_rect(color="black")) +
  labs(y = TeX("$X_2$"),
       x = TeX("$X_1$"),
       pch = TeX(" $Z$"))

ggsave("figures/show_cem_calipers_nocsm.png", width=6, height=3.5)




