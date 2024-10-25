

require(mvtnorm)
library( tidyverse )

# build gradient-looking background ---------------------------------------
MU <- c(1.5,1.5)
SIG <- matrix(c(1,0.5,0.5,1), nrow=2)

toy_background <- expand.grid(
  x1 = seq(0.025,3,by=0.01),
  x2 = seq(0.025,3,by=0.01)
) %>%
  rowwise() %>%
  mutate(z = dmvnorm(x = c(x1,x2),
                     mean = MU,
                     sigma = SIG))

saveRDS(toy_background,
        "data/toy_background_df.rds")
