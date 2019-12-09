
#'@title Estimate CAUSE Nuisance Parameters
#'@description Estimates bivariate distribution of summary statistics and rho, the
#'correlation between summary statistics due to overlapping samples or population structure.
#'@param X An object of class cause_data containing data for the two traits.
#'@param variants A vector of variants to include. This list should be approximately LD
#'pruned and include variants genome wide.
#'@return An object of class cause_params
#'@export
est_cause_params <- function(X, variants, optmethod = c("mixSQP", "mixIP"),
                             null_wt = 10, max_candidates=Inf){
  optmethod <- match.arg(optmethod)
  stopifnot(inherits(X, "cause_data"))
  if(!all(variants %in% X$snp)){
    warning("Warning: Not all `variants` are in data.", call.=FALSE)
  }
  X <- filter(X, snp %in% variants)
  X <- new_cause_data(X)
  if(nrow(X) < 1e5){
    warning("Fewer than 100,000 variants are being used to estimate parametrs. ",
            "This can cause problems and is not recomended. You are using ", nrow(X),
            " variants.\n")
  }
  cat("Estimating CAUSE parameters with ", nrow(X), " variants.\n")
  mix_grid <- variance_pair_candidates(X, optmethod=optmethod, max_candidates=max_candidates)

  params <- map_pi_rho(X, mix_grid, optmethod=optmethod, null_wt = null_wt)
  #Filter out grid points with low mixing proportion
  params$mix_grid <- dplyr::filter(params$mix_grid, zapsmall(pi) > 0)

  return(params)
}


#'@export
variance_pair_candidates <- function(X, optmethod = c("mixSQP", "mixIP",
                                                      "cxxMixSquarem", "mixEM",
                                                      "mixVBEM", "w_mixEM"),
                                     gridmult=sqrt(2), max_candidates = Inf){
  stopifnot(inherits(X, "cause_data"))
  optmethod <- match.arg(optmethod)

  get_candidates1 <- function(X, optmethod, gridmult){
    fit1 <- with(X, ash(betahat = beta_hat_1, sebetahat = seb1,
                        mixcompdist = "normal", prior="nullbiased",
                        optmethod=optmethod, gridmult=gridmult))
    with(fit1$fitted_g, sd[!zapsmall(pi)==0])
  }

  get_candidates2 <- function(X, optmethod, gridmult){
    fit2 <- with(X, ash(betahat = beta_hat_2, sebetahat = seb2,
                      mixcompdist = "normal", prior="nullbiased",
                      optmethod = optmethod, gridmult=gridmult))
    with(fit2$fitted_g, sd[!zapsmall(pi)==0])
  }
  gridmult1 <- gridmult2 <- gridmult

  sigma1 <- get_candidates1(X, optmethod, gridmult)

  while(length(sigma1) > max_candidates-2){
    gridmult1 <- sqrt((gridmult1^2)*2)
    #cat(gridmult1^2, "  ")
    sigma1 <- get_candidates1(X, optmethod, gridmult1)
    #cat(length(sigma1), "\n")
  }

  sigma2 <- get_candidates2(X, optmethod, gridmult)
  while(length(sigma2) > max_candidates-2){
    gridmult2 <- sqrt((gridmult2^2)*2)
    #cat(gridmult2^2, "  ")
    sigma2 <- get_candidates2(X, optmethod, gridmult2)
    #cat(length(sigma2), "\n")
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


