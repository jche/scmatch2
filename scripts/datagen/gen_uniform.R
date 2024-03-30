gen_uniform<- function(n = 500,
                       p_trt = 0.5,
                      xlim = c(0, 2.25),
                      ylim = c(0, 2.25),
                      MU = c(1.5,1.5),
                      SIG = matrix(c(1,0.5,0.5,1), nrow=2)) {
  uniform_df<- data.frame(
    X1=runif(n, min=xlim[1], max=xlim[2]),
    X2=runif(n, min=ylim[1], max=ylim[2]),
    Z = sample(c(0,1), size=n,
           replace = T,
           prob=c(1-p_trt, p_trt))
  )


  uniform_df <-
    uniform_df %>%
    rowwise() %>%
    mutate(Y0 =  mvtnorm::dmvnorm(x = c(X1,X2),
                                            mean = MU,
                                            sigma = SIG)  )
  return(uniform_df)
}
