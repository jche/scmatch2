
# functions for simulating data

# require(aciccomp2016)
library(mvtnorm)

# toy example -------------------------------------------------------------
#' Generate covariates X1, X2 for the toy example
#'
#' @param n Number of samples
#' @param X1_ctrs A vector of two elements, representing the two different
#'  X1 locations
#' @param X2_ctrs Same as X1_ctrs, but
#' @param SD Std. deviation of the
#'
#' @return A tibble object with two columns X1, X2
#' @export
#'
#' @examples n = 100
#' X1_ctrs <- c(0.25, 0.75)
#' X2_ctrs <- c(0.25, 0.75)
#' SD <- 0.1
#' res <- gen_toy_covar(n, X1_ctrs, X2_ctrs, SD)
gen_toy_covar <- function(n, X1_ctrs, X2_ctrs, SD){
  tibble::tibble(
    X1 = c(rnorm(n/2, mean=X1_ctrs[1], sd=SD),
           rnorm(n/2, mean=X1_ctrs[2], sd=SD)),
    X2 = c(rnorm(n/2, mean=X2_ctrs[1], sd=SD),
           rnorm(n/2, mean=X2_ctrs[2], sd=SD))
  )
}

# generate data mostly in grid (0,1) x (0,1)
#' Generate toy data
#'
#' @param nc (something)
#' @param nt (something)
#' @param f0_sd (something)
#' @param f0_fun (something)
#' @param tx_effect_fun (something)
#' @param ctr_dist (something)
#' @param prop_nc_unif (something)
#'
#' @return A tibble containing a toy dataset
#' @export
#'
gen_df_adv <- function(nc, nt,
                       f0_sd = 0.1,   # homoskedastic noise
                       f0_fun = function(X1, X2) {1},
                       tx_effect_fun = function(X1, X2) {1},
                       ctr_dist = 0.5, # from 0 to 1;
                       #            lower means better overlap
                       prop_nc_unif = 1/3 # proportion of uniform controls
                       #            lower means worse overlap
) {
  require( tidyverse )

  SD <- 0.1

  c1 <- 0.5 - ctr_dist/2
  c2 <- 0.5 + ctr_dist/2

  # nt tx units clustered at (0.25,0.25) and (0.75,0.75)
  #     which translates to X1_ctrs=(0.25,0.75) and X2_ctrs=(0.25,0.75)
  dat_txblobs <-
    gen_toy_covar(nt, X1_ctrs=c(c1,c2), X2_ctrs=c(c1,c2), SD)
  dat_txblobs$Z <- T

  nc_unif <- ceiling(nc * prop_nc_unif)
  # (nc-nc_unif) co units clustered at (0.75,0.25) and (0.25,0.75)
  #     which translates to X1_ctrs=(0.75,0.25) and X2_ctrs=(0.25,0.75)
  dat_coblobs <-
    gen_toy_covar(
      nc-nc_unif,
      X1_ctrs=c(c2,c1),
      X2_ctrs=c(c1,c2),
      SD)
  dat_coblobs$Z <- F

  # nc_unif co units uniformly scattered on (0,1) box to give some
  #   randomly good controls
  dat_conear <- tibble(
    X1 = runif(nc_unif),
    X2 = runif(nc_unif),
    Z  = F
  )

  dat <- bind_rows(dat_txblobs,
                   dat_coblobs,
                   dat_conear)

  res <- dat %>%
    mutate(noise = rnorm(n(), mean=0, sd=f0_sd)) %>%
    mutate(Y0 = f0_fun(X1,X2) + noise) %>%
    mutate(Y1 = Y0 + tx_effect_fun(X1,X2),
           Y  = ifelse(Z, Y1, Y0)) %>%
    mutate(id = 1:n(), .before=X1)
  # print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))

  return(res)
}

if (F) {
  df <- gen_df_adv(
    nc=nc,
    nt=nt,
    f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) {1},
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>%
        mvtnorm::dmvnorm(mean = c(0.5,0.5),
                         sigma = matrix(c(1,0.8,0.8,1), nrow=2))}) %>%
    mutate(X3 = X1*X2)
  df %>%
    ggplot(aes(X1,X2)) +
    geom_point(aes(color=Y0)) +
    facet_wrap(~Z)
}


#' Generate k-dimensional blobs around specified centers
#' (Helper function replacing the original 2D gen_toy_covar)
#'
#' @param n Number of points to generate.
#' @param centers A list of k-dimensional numeric vectors representing cluster centers.
#' @param k Dimensionality.
#' @param sd Standard deviation for the blobs (assumed spherical).
#' @return A tibble with k columns (X1, ..., Xk).
#'
gen_toy_covar_k <- function(n, centers, k, sd) {
  n_centers <- length(centers)
  # Return empty tibble with correct structure if n=0 or no centers
  if (n <= 0 || n_centers == 0) {
    col_names <- paste0("X", 1:k)
    named_list <- setNames(replicate(k, numeric(), simplify = FALSE), col_names)
    return(tibble(!!!named_list))
  }

  # Check center dimensions
  if (!all(sapply(centers, length) == k)) {
    stop("All centers must have dimension k.")
  }

  # Assign points to centers roughly evenly
  assignments <- sample(1:n_centers, size = n, replace = TRUE)
  counts <- table(factor(assignments, levels = 1:n_centers))

  # Generate points for each center using multivariate normal
  all_points_list <- vector("list", n_centers)
  cov_matrix <- diag(sd^2, k) # Simple spherical covariance

  for (i in 1:n_centers) {
    count_i <- as.integer(counts[i]) # Get count for center i
    if (!is.na(count_i) && count_i > 0) {
      center_vec <- centers[[i]]
      points <- mvtnorm::rmvnorm(count_i, mean = center_vec, sigma = cov_matrix)
      # Ensure points is a matrix even if count_i is 1
      if(count_i == 1) points <- matrix(points, nrow=1)

      # Add column names to the matrix before converting to tibble
      colnames(points) <- paste0("X", 1:k)

      # Now convert to tibble with proper names
      all_points_list[[i]] <- as_tibble(points)
    }
  }

  # Then remove the rename operation later since columns are already named
  df <- bind_rows(all_points_list)

  return(df)
}


#' Generate k-dimensional toy data (generalized version)
#'
#' @param nc Number of control units.
#' @param nt Number of treated units.
#' @param k Dimensionality of the feature space.
#' @param f0_sd Standard deviation of the noise term.
#' @param f0_fun Function for baseline potential outcome Y0 (expects k-dim matrix X).
#' @param tx_effect_fun Function for treatment effect (expects k-dim matrix X).
#' @param ctr_dist Distance parameter influencing cluster separation.
#' @param prop_nc_unif Proportion of control units drawn uniformly from [0,1]^k box.
#'
#' @return A tibble containing a k-dimensional toy dataset.
#' @export
#'
gen_df_adv_k <- function(nc, nt, k, # Added k
                         f0_sd = 0.1,
                         # Default functions now expect matrix X
                         f0_fun = function(X) { rep(1, nrow(X)) },
                         tx_effect_fun = function(X) { rep(1, nrow(X)) },
                         ctr_dist = 0.5,
                         prop_nc_unif = 1/3
) {

  SD <- 0.1 # Base SD for blobs

  # Define k-dimensional centers
  c1 <- 0.5 - ctr_dist / 2
  c2 <- 0.5 + ctr_dist / 2

  # Tx centers: e.g., (c1, c1, ...) and (c2, c2, ...)
  tx_centers <- list(rep(c1, k), rep(c2, k))

  # Control centers: e.g., swap first dim: (c2, c1, ...) and (c1, c2, ...)
  co_centers <- list(c(c2, rep(c1, max(0, k - 1))), c(c1, rep(c2, max(0, k - 1))))
  # Handle k=1 case explicitly if needed (though max handles it)
  # if (k == 1) { co_centers <- list(c(c2), c(c1)) }

  # Generate treated blobs
  dat_txblobs <- gen_toy_covar_k(nt, centers = tx_centers, k = k, sd = SD)
  if (nrow(dat_txblobs) > 0) dat_txblobs$Z <- TRUE

  # Generate control blobs
  nc_blob <- nc - ceiling(nc * prop_nc_unif)
  dat_coblobs <- gen_toy_covar_k(nc_blob, centers = co_centers, k = k, sd = SD)
  if (nrow(dat_coblobs) > 0) dat_coblobs$Z <- FALSE

  # Generate uniform controls in [0,1]^k
  nc_unif <- ceiling(nc * prop_nc_unif)
  if (nc_unif > 0) {
    unif_matrix <- matrix(runif(nc_unif * k), ncol = k)
    dat_conear <- as_tibble(unif_matrix)
    colnames(dat_conear) <- paste0("X", 1:k)
    dat_conear$Z <- FALSE
  } else {
    # Create an empty tibble with correct columns if nc_unif is 0
    col_names <- paste0("X", 1:k)
    named_list <- setNames(replicate(k, numeric(), simplify = FALSE), col_names)
    dat_conear <- tibble(!!!named_list, Z = logical())
  }

  # Combine data
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)

  # --- Crucial change for function calling ---
  # Select the X columns into a matrix to pass to functions
  X_matrix <- dat %>% select(all_of(paste0("X", 1:k))) %>% as.matrix()
  # ---

  # Apply functions using the prepared X_matrix
  res <- dat %>%
    mutate(noise = rnorm(n(), mean = 0, sd = f0_sd)) %>%
    # Pass the matrix X_matrix to the functions
    mutate(
      Y0_denoised = f0_fun(X_matrix),
      Y0 = f0_fun(X_matrix) + noise
    ) %>%
    mutate(
      Y1_denoised = Y0_denoised + tx_effect_fun(X_matrix),
      Y1 = Y0 + tx_effect_fun(X_matrix)
    ) %>%
    mutate(
      Y = ifelse(Z, Y1, Y0),
      Y_denoised = ifelse(Z, Y1_denoised, Y0_denoised)
    ) %>%
    # Ensure id is added correctly (use first X col name if k>=1)
    mutate(id = 1:n(), .before = 1)

  return(res)
}


#' Generate a toy dataset with a single treatment effect (Updated to use gen_df_adv_k)
#'
#' @param k Dimensionality of the feature space. Default is 2.
#' @param nc Number of control units. Default is 500.
#' @param nt Number of treated units. Default is 100.
#' @param ctr_dist Distance between the two control clusters. Default is 0.5.
#' @param prop_nc_unif Proportion of control units drawn from a uniform distribution. Default is 1/3.
#' @param f0_sd Standard deviation for the noise in the potential outcome function f0. Default is 0.5.
#'
#' @return A tibble with the toy dataset
#' @export
#'
gen_one_toy <- function( k = 2, # Added dimensionality parameter k, default 2
                         nc = 500, nt = 100,
                         ctr_dist = 0.5,
                         prop_nc_unif = 1/3,
                         f0_sd = 0.5 ){

  # --- Input validation ---
  if ( nc < 5 ) {
    warning( "Very small control group in gen_one_toy!" )
  }
  if (!is.numeric(k) || k < 1 || floor(k) != k) {
    stop("k must be a positive integer.")
  }
  # --- End Input validation ---

  # --- Define k-dimensional functions ---
  # These functions expect a matrix 'X'
  tx_effect_fun <- function(X) {
    if (!is.matrix(X) && !is.data.frame(X)) stop("tx_effect_fun expects a matrix/df.")
    if (ncol(X) < k) stop(paste("tx_effect_fun requires input with at least", k, "columns."))
    rowSums(X[, 1:k, drop = FALSE] * 3)
  }

  f0_fun <- function(X) {
    if (!is.matrix(X) && !is.data.frame(X)) stop("f0_fun expects a matrix/df.")
    if (ncol(X) < k) stop(paste("f0_fun requires input with at least", k, "columns."))
    mean_vec <- rep(0.5, k)
    sigma_mat <- matrix(0.8, nrow = k, ncol = k)
    diag(sigma_mat) <- 1
    densities <- mvtnorm::dmvnorm(X[, 1:k, drop=FALSE], mean = mean_vec, sigma = sigma_mat)
    densities * 20
  }
  # --- End function definitions ---

  # --- Call the NEW k-dimensional generator function ---
  # Pass k explicitly now
  gen_df_adv_k(
    nc = nc,
    nt = nt,
    k = k, # Pass k
    f0_sd = f0_sd,
    tx_effect_fun = tx_effect_fun,
    f0_fun = f0_fun,
    ctr_dist = ctr_dist,
    prop_nc_unif = prop_nc_unif
  )
  # --- End function call ---
}


#'
#' #' Generate a toy dataset with a single treatment effect
#' #'
#' #' @param ctr_dist Distance between the two control clusters
#' #'
#' #' @return A tibble with the toy dataset
#' #' @export
#' #'
#' gen_one_toy <- function( nc = 500, nt = 100,
#'                          ctr_dist = 0.5,
#'                          prop_nc_unif = 1/3,
#'                          f0_sd = 0.5 ){
#'   if ( nc < 5 ) {
#'     warning( "Very small control group in gen_one_toy!" )
#'   }
#'   gen_df_adv(
#'     nc=nc,
#'     nt=nt,
#'     f0_sd = f0_sd,
#'     tx_effect_fun = function(X1, X2) {3*X1+3*X2},
#'     f0_fun = function(x,y) {
#'       matrix(c(x,y), ncol=2) %>%
#'         mvtnorm::dmvnorm( mean = c(0.5,0.5),
#'                           sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
#'     },
#'     ctr_dist = ctr_dist,
#'     prop_nc_unif = prop_nc_unif
#'   )
#' }



# canonical examples ------------------------------------------------------

# generate sample dataset from Hainmueller (2012), exactly
#  - note: generates a big population and samples nt/nco units
# paper settings:
#  - n (total number of units): 300, 600, 1500
#  - r (ratio of co/tx units): 1,2,5
#  - sigma_e: n30 = low overlap, n100 = high overlap, chi5 = weird overlap
#  - outcome: "linear", "nl1", "nl2"
#  - sigma_y: 1
gen_df_hain <- function(nt = 50,
                        nc = 250,
                        sigma_e = c("chi5", "n30", "n100"),
                        outcome = c("linear", "nl1", "nl2"),
                        sigma_y = 1,
                        ATE = 0) {
  sigma_e <- match.arg(sigma_e)
  outcome <- match.arg(outcome)
  NUMSAMP <- max(nt, nc)*10

  if (sigma_e == "n30") {
    eps_e <- rnorm(NUMSAMP, sd=sqrt(30))
  } else if (sigma_e == "n100"){
    eps_e <- rnorm(NUMSAMP, sd=sqrt(100))
  } else if (sigma_e == "chi5") {
    eps_e <- rchisq(NUMSAMP, df=5)
    eps_e <- (eps_e - mean(eps_e)) / sd(eps_e)
    eps_e <- eps_e * sqrt(67.6) + 0.5
  }

  df <- as_tibble(rmvnorm(NUMSAMP,
                          mean = c(0,0,0),
                          sigma = matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1), ncol=3))) %>%
    mutate(V4 = runif(NUMSAMP, min=-3, max=3),
           V5 = rchisq(NUMSAMP, df=1),
           V6 = sample(c(T,F), NUMSAMP, replace=T, prob=c(0.5,0.5)),
           Z  = (V1 + 2*V2 - 2*V3 + V4 - 0.5*V5 + V6 + eps_e) > 0,
           id = 1:n()) %>%
    relocate(id)

  if (outcome == "linear") {
    df <- df %>%
      mutate(Y0 = V1 + V2 + V3 - V4 + V5 + V6 + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl1") {
    df <- df %>%
      mutate(Y0 = V1 + V2 + 0.2*V3*V4 - sqrt(V5) + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl2") {
    df <- df %>%
      mutate(Y0 = (V1 + V2 + V5)^2 + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl3") {   # not in paper
    df <- df %>%
      mutate(Y0 = V1*V2 + V2*V5 + V1*V5 + 3*cos(V3) + 2*sin(V2) +
               V3*V4 + 0.5*V4^2 + rnorm(NUMSAMP, sd=sigma_y))
  }
  df <- df %>%
    mutate(Y1 = Y0 + ATE,
           Y  = ifelse(Z, Y1, Y0))

  # get correct number of observations
  if (nt > sum(df$Z) | nc > nrow(df)-sum(df$Z)) {
    warning("Insufficient population pool: increase NUMSAMP")
  }
  df_tx <- df %>%
    filter(Z) %>%
    slice_sample(n=nt)
  df_co <- df %>%
    filter(!Z) %>%
    slice_sample(n=nc)

  # bind rows, tx then co
  df_final <- df_tx %>%
    bind_rows(df_co)
  names(df_final) <- names(df_final) %>%
    str_replace("V", "X")

  return(df_final)
}

gen_df_acic <- function(model.trt="step",
                        root.trt=0.35,
                        overlap.trt="full",
                        model.rsp="linear",
                        alignment=0.75,
                        te.hetero="high",
                        random.seed=1,
                        n=1000, p=10) {
  if (!requireNamespace("aciccomp2016", quietly = TRUE)) {
    stop("Package 'aciccomp2016' is needed for this function. Please install it with: install.packages('kbal')",
         call. = FALSE)
  }
  # idea: input_2016 is a 4802 x 58 df of COVARIATES
  #  - dgp_2016() outputs a dataset of POs for input_2016
  #  - see parameters_2016 for example inputs

  # rows:
  #  - model.trt: tx assignment model (linear, polynomial, step)
  #  - root.trt: percentage of tx units (any #; 0.35, 0.65)
  #  - overlap.trt: overlap between tx/co (one-term, full)
  #  - model.rsp: outcome model (linear, exponential, step)
  #  - alignment: prop of tx assignment covs used in response fxn (any #; 0, 0.25, 0.75)
  #  - te.hetero: tx-effect heterogeneity (none, med, high)

  # NOTE: for some god-forsaken reason, subsetting rows breaks things
  #  - i.e., it makes the browser() functionality start running...
  #  - the browser() call is from newDiscreteBaseFunction(),
  #    where numLevels <= 1...?

  set.seed(random.seed)

  # only subset numeric columns
  numeric_cols <- input_2016 %>%
    summarize(across(everything(), ~length(unique(.)))) %>%
    pivot_longer(everything()) %>%
    filter(value >= 20) %>%
    pull(name)
  base_df <- input_2016[,numeric_cols]

  if (n > nrow(base_df)) {
    error(glue("n too large; needs to be <= {nrow(base_df)}"))
  }
  if (p > ncol(base_df)) {
    error(glue("p too large; needs to be <= {ncol(base_df)}"))
  }
  base_df <- base_df[,1:p] %>%
    sample_n(n, replace=F)

  # # NECESSARY: ensure that there aren't any factor levels without entries
  # base_df$x_2 <- as.factor(as.character(base_df$x_2))

  # base_df %>%
  #   summarize(across(everything(), ~length(unique(.))))
  # browser()
  df <-
    dgp_2016(base_df,
             parameters = list(
               model.trt = model.trt,
               root.trt = root.trt,
               overlap.trt = overlap.trt,
               model.rsp = model.rsp,
               alignment = alignment,
               te.hetero = te.hetero),
             random.seed = random.seed)

  df <- bind_cols(base_df, df) %>%
    rename_with(.cols = everything(), toupper) %>%
    rename_with(.cols = everything(), ~str_remove(., "_")) %>%
    rename_with(.cols = everything(), ~str_remove(., "\\.")) %>%
    mutate(id = 1:n(), .before=X1,
           Z = Z==1) %>%
    as_tibble()

  return(df)
}

gen_df_kang <- function(n=1000) {
  as_tibble(rmvnorm(n, mean=rep(0,4), sigma=diag(4))) %>%
    mutate(Y = 210 + 27.4*V1 + 13.7*(V2+V3+V4) + rnorm(n),
           e = rje::expit(-V1 + 0.5*V2 - 0.25*V3 - 0.1*V4)) %>%
    mutate(X1 = exp(V1/2),
           X2 = V2/(1+exp(V1)) + 10,
           X3 = (V1*V3/25 + 0.6)^3,
           X4 = (V2 + V4 + 20)^2) %>%
    rowwise() %>%
    mutate(Z = sample(c(T,F), size=1, replace=T,
                      prob=c(e,1-e))) %>%
    ungroup() %>%
    mutate(id = 1:n, .before=V1)
}


if (F) {
  gen_df_acic(n=1000, p=10)

  input_2016[,1:10] %>%
    sample_n(1000,replace=F) %>%
    dgp_2016(
      parameters = list(
        model.trt="linear",
        root.trt=0.35,
        overlap.trt="full",
        model.rsp="linear",
        alignment=0.75,
        te.hetero="none"),
      random.seed=1)
}




