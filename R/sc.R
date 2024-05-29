
# functions related to synthetic controls


#' Solve the synth QP directly (code taken from augsynth package)
#'
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


# n0 + 1 decision variables: weight for each co unit & one slack variable
synth_lp <- function(X1, X0, V) {

  n0 <- nrow(X0)
  p  <- ncol(X0)

  # build (p x n0) constant matrix
  VX0T <- V %*% t(X0)
  con_mat <- rbind(
    c(0, rep(1, n0)),
    cbind(rep(-1, p), -VX0T),
    cbind(rep(-1, p), VX0T)
  )
  # build (p x 1) constraint vector
  VX1T <- V %*% X1

  obj <- c(1, rep(0,n0))
  dir <- c("=", rep("<=", 2*p))
  rhs <- c(1, -VX1T, VX1T)

  # browser()
  sol <- lpSolve::lp(objective.in = obj,
                     const.mat = con_mat,
                     const.dir = dir,
                     const.rhs = rhs)

  return(sol$solution[-1])
}


# main function: generates SC weights
gen_sc_weights <- function(d,
                           match_cols,
                           dist_scaling,
                           metric = c("maximum", "euclidean", "manhattan")) {

  metric <- match.arg(metric)

  # edge cases
  if (nrow(d) == 0) {
    return(tibble())
  } else if (nrow(d) == 2) {
    return(d %>%
             mutate(unit = c("tx1", "c1"),
                    weights = c(1,1)))
  }
  X1 = d[1,] %>%
    select(all_of(match_cols)) %>%
    as.numeric()
  X0 = d[-1,] %>%
    select(all_of(match_cols)) %>%
    as.matrix()
  if (metric == "maximum") {
    # run linear program
    if ( is.null(nrow(dist_scaling)) ||  nrow(dist_scaling==1) ){
        V = diag(dist_scaling, length(match_cols))
    }else{
      V = dist_scaling %>%
        select(all_of(match_cols)) %>%
        as.numeric() %>%
        diag()
    }


    sol <- synth_lp(X1 = X1,
                    X0 =X0,
                    V  = V)
  } else if (metric == "euclidean") {
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
  } else if (metric == "manhattan") {
    stop("Linear program for L1-distance minimization is not currently implemented.")
  }

  d %>%
    mutate(unit = c("tx1", paste0("c", 1:(n()-1))),
           weights = c(1, sol))
}






