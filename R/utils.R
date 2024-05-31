
# helpful utility functions

logit <- function(x) {
  log(x/(1-x))
}
invlogit <- function(x) {
  exp(x) / (1+exp(x))
}

expit <- function(x) {
  1 / (1+exp(-x))
}

rmse <- function(resid){
  sqrt(sum(resid^2))
}

weighted_var <- function(x, wt) {
  n <- length(x)
  wm <- weighted.mean(x, wt)
  sum(wt * (x - wm)^2) * n / (n-1)
}

weighted_se <- function(x, wt) {
  sqrt(weighted_var(x, wt) / length(x))
}

get_x_vars <- function(df) {
  names(df) %>%
    grep("^X", ., value = TRUE)
}
