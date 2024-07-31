library(tidyverse)

gen_six_points <- function(x1 = 0.75, x2 = 2.25,
                           y1 = 0.75, y2 = 2.25,
                           NUDGE = 0.125,
                           tau=function(x1,x2) 3*x1+3*x2,
                           scale_Y0 = 1){
  tx_dat <- tibble(
    X1 = c(x1, x2),
    X2 = c(y1, y2),
    Z  = 1,
  )
  co_dat <- tibble(
    X1 = c(x1, x2),
    X2 = c(y2, y1),
    Z  = 0
  )
  four_points_df <- rbind(tx_dat, co_dat)

  local_controls <- tibble(
    X1 = c(x1-1.5*NUDGE, x2-1.5*NUDGE),
    X2 = c(x1-0.5*NUDGE, x2+1.5*NUDGE),
    Z  = 0
  )

  six_points_df <-
    rbind(four_points_df, local_controls)


  MU = c(1.5,1.5)
  SIG = matrix(c(1,0.5,0.5,1), nrow=2)

  six_points_df <-
    six_points_df %>%
    rowwise() %>%
    mutate(Y0 = scale_Y0 * mvtnorm::dmvnorm(x = c(X1,X2),
                       mean = MU,
                       sigma = SIG)  ) %>%
    mutate(Y1 = Y0 + tau(X1,X2) ) %>%
    mutate(Y = (1-Z)*Y0 + Z*Y1)
  return(six_points_df)
}

if ( FALSE ) {
  ggplot( gen_six_points(), aes( X1, X2 ) ) +
    geom_point(aes(color=Z), size=5)
}
