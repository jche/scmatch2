
# plot_toy.R
# Collection of plotting functions


library(ggplot2)
require(latex2exp)

create_background_plot <- function(background_df){
  background_plot <-
    ggplot(background_df,aes(x1,x2)) +
    geom_tile(aes(fill=z),
              alpha=0.65)+
    scale_fill_gradient(low="blue", high="orange") +
    geom_contour(aes(z=z),
                 color="black",
                 alpha = 0.1,
                 lty = "solid",
                 bins = 5)
  return(background_plot)
}

create_toy_example_four_points_plot <-
  function(background_plot, four_points_df){
    points_on_background_plot <-
      background_plot+
      geom_point(data=four_points_df,
                 aes(pch=as.factor(z)),
                 size=2)

    with_label_plot <-
      points_on_background_plot +
      annotate("text", x=X1+NUDGE, y=Y1+NUDGE, label="t[1]", parse=T) +
      annotate("text", x=X2+NUDGE, y=Y2+NUDGE, label="t[2]", parse=T) +
      annotate("text", x=X1+NUDGE, y=Y2+NUDGE, label="c[1]", parse=T) +
      annotate("text", x=X2+NUDGE, y=Y1+NUDGE, label="c[2]", parse=T)

    with_axis_lines_plot <-
      with_label_plot +
      geom_segment(aes(x=0, y=0, xend=0, yend=3),
                   arrow = arrow(length=unit(0.2, "cm")),
                   linewidth=1) +
      geom_segment(aes(x=0, y=0, xend=3, yend=0),
                   arrow = arrow(length=unit(0.2, "cm")),
                   linewidth=1)

    with_theme_plot <-
      with_axis_lines_plot +
      theme(
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank()
      ) +
      labs(y = TeX("$X_2$"),
           x = TeX("$X_1$"),
           pch = TeX(" $Z$"))

    shifted_plot <-
      with_theme_plot +
      coord_cartesian(xlim = c(0.1,3),
                      ylim = c(0.05,3))

    return(shifted_plot)
}

create_toy_example_six_points_plot <-
  function(background_plot, six_points_df, draw_circle=T){
    toy_example_four_points_plot <-
      create_toy_example_four_points_plot(background_plot, six_points_df)
    with_point_text_plot <-
      toy_example_four_points_plot +
      annotate("text",
               x=X1-2.9*NUDGE,
               y=Y1-0.5*NUDGE,
               label="c[3]",
               parse=T) +
      annotate("text",
               x=X2-2.3*NUDGE,
               y=Y2+0.7*NUDGE,
               label="c[4]",
               parse=T)
    if (draw_circle){
      with_caliper_plot <-
        with_point_text_plot +
        ggforce::geom_circle(aes(x0=X1, y0=Y1, r=NUDGE*4),
                             lty = "dotted") +
        ggforce::geom_circle(aes(x0=X2, y0=Y2, r=NUDGE*4),
                             lty = "dotted")
    }else{
      with_caliper_plot <-
        with_point_text_plot
    }

    return(with_caliper_plot)
}
