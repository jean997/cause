
#'@title Estimate CAUSE Nuisance Parameters
#'@description Estimates bivariate distribution of summary statistics and rho, the
#'correlation between summary statistics due to overlapping samples or population structure.
#'@param X An object of class cause_data containing data for the two traits.
#'@param variants A vector of variants to include. This list should be approximately LD
#'pruned and include variants genome wide.
#'@return An object of class cause_params
#'@export
est_cause_params <- function(X, variants, optmethod = c("mixSQP", "mixIP"),
                             null_wt = 10, max_candidates=Inf, control=list()){
  optmethod <- match.arg(optmethod)
  if(optmethod == "mixSQP"){
     control0 <- list(verbose = FALSE, eps = 1e-08, numiter.em = 10, maxiter.activeset = 100, tol.svd = 0)
     control <- modifyList(control0, control, keep.null = TRUE)
  }
  if(!inherits(X, "cause_data")){
    X <- new_cause_data(X)
  }
  if(!all(variants %in% X$snp)){
    warning("Warning: Not all `variants` are in data.", call.=FALSE)
  }
  X <- filter(X, snp %in% variants) %>%
       new_cause_data()
  if(nrow(X) < 1e5){
    warning("Fewer than 100,000 variants are being used to estimate parametrs. ",
            "This can cause problems and is not recomended. You are using ", nrow(X),
            " variants.\n")
  }
  cat("Estimating CAUSE parameters with ", nrow(X), " variants.\n")
  mix_grid <- variance_pair_candidates(X, optmethod=optmethod, max_candidates=max_candidates, control=control)

  params <- map_pi_rho(X, mix_grid, optmethod=optmethod, null_wt = null_wt, control=control, warm=TRUE)
  #Filter out grid points with low mixing proportion
  params$mix_grid <- dplyr::filter(params$mix_grid, zapsmall(pi) > 0)

  return(params)
}


#'@export
variance_pair_candidates <- function(X, optmethod = c("mixSQP", "mixIP"),
                                     gridmult=sqrt(2), max_candidates = Inf, control=list()){
  stopifnot(inherits(X, "cause_data"))
  optmethod <- match.arg(optmethod)

  get_candidates <- function(beta_hat, se_beta_hat, optmethod, gridmult, control){
    fit <- ash(betahat = beta_hat, sebetahat = se_beta_hat,
                        mixcompdist = "normal", prior="nullbiased",
                        optmethod=optmethod, gridmult=gridmult, control=control)
    with(fit$fitted_g, sd[!zapsmall(pi)==0])
  }

  gridmult1 <- gridmult2 <- gridmult

  sigma1 <- get_candidates(X$beta_hat_1, X$seb1, optmethod, gridmult, control)
  while(length(sigma1) > max_candidates-2){
    gridmult1 <- sqrt((gridmult1^2)*2)
    sigma1 <- get_candidates(X$beta_hat_1, X$seb1, optmethod, gridmult1, control)
  }

  sigma2 <- get_candidates(X$beta_hat_2, X$seb2, optmethod, gridmult, control)
  while(length(sigma2) > max_candidates-2){
    gridmult2 <- sqrt((gridmult2^2)*2)
    sigma2 <- get_candidates(X$beta_hat_2, X$seb2, optmethod, gridmult2, control)
  }

  sigma1 <- c(sigma1, 2*max(sigma1), 4*max(sigma1))
  sigma2 <- c(sigma2, 2*max(sigma2), 4*max(sigma2))
  mix_grid <- expand.grid("S1"=sigma1, "S2"=sigma2, "pi"=0)
  new_cause_grid(mix_grid)
}

new_cause_grid <- function(x = data.frame()){
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c("S1", "S2", "pi") %in% names(x)))
  structure(x, class = c("cause_grid", "data.frame"))
}


