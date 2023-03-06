
# functions related to synthetic controls


#' Solve the synth QP directly (from augsynth package)
#' @param X1 Target vector
#' @param X0 Matrix of control outcomes
#' @param V Scaling matrix
synth_qp <- function(X1, X0, V) {
  
  Pmat <- X0 %*% V %*% t(X0)
  qvec <- - t(X1) %*% V %*% t(X0)
  
  n0 <- nrow(X0)
  A <- rbind(rep(1, n0), diag(n0))
  l <- c(1, numeric(n0))
  u <- c(1, rep(1, n0))
  
  # browser()
  
  settings = osqp::osqpSettings(verbose = FALSE,
                                eps_rel = 1e-8,
                                eps_abs = 1e-8)
  sol <- osqp::solve_osqp(P = Pmat, q = qvec,
                          A = A, l = l, u = u, 
                          pars = settings)
  
  return(sol$x)
}


# main function: generates SC weights
gen_sc_weights <- function(d, match_cols, dist_scaling) {
  
  # edge cases
  if (nrow(d) == 0) {
    return(tibble())
  } else if (nrow(d) == 2) {
    return(d %>% 
             mutate(unit = c("tx1", "c1"),
                    weights = c(1,1)))
  }
  # browser()
  
  # run osqp
  #  - note: square V, since we use (V^T V) within the euclidean distance!
  sol <- synth_qp(X1 = d[1,] %>% 
                    select(all_of(match_cols)) %>% 
                    as.numeric(),
                  X0 = d[-1,] %>% 
                    select(all_of(match_cols)) %>% 
                    as.matrix(),
                  V  = dist_scaling %>% 
                    select(all_of(match_cols)) %>% 
                    as.numeric() %>% 
                    map_dbl(~.x^2) %>% 
                    diag())
  
  d %>% 
    mutate(unit = c("tx1", paste0("c", 2:(n()))),
           weights = c(1, sol))
}






