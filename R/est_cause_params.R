
#'@title Estimate CAUSE Nuisance Parameters
#'@description Estimates bivariate distribution of summary statistics and rho, the
#'correlation between summary statistics due to overlapping samples or population structure.
#'@param X An object of class cause_data containing data for the two traits.
#'@param variants A vector of variants to include. This list should be approximately LD
#'pruned and include variants genome wide.
#'@return An object of class cause_params
#'@export
est_cause_params <- function(X, variants){
  stopifnot(inherits(X, "cause_data"))
  if(!all(variants %in% X$snp)){
    stop("Not all `variants` are in data.", call.=FALSE)
  }
  X <- filter(X, snp %in% variants)
  X <- new_cause_data(X)
  if(nrow(X) < 1e5){
    warning("Fewer than 100,000 variants are being used to estimate parametrs. ",
            "This can cause problems and is not recomended. You are using ", nrow(X),
            " variants.\n")
  }
  mix_grid <- variance_pair_candidates(X)

  params <- map_pi_rho(X, mix_grid)

  return(params)
}


#'@export
variance_pair_candidates <- function(X, optmethod = c("mixIP", "mixSQP",
                                                      "cxxMixSquarem", "mixEM",
                                                      "mixVBEM", "w_mixEM")){
  stopifnot(inherits(X, "cause_data"))
  optmethod <- match.arg(optmethod)
  fit1 <- with(X, ash(betahat = beta_hat_1, sebetahat = seb1,
                       mixcompdist = "normal", prior="nullbiased",
                      optmethod=optmethod))
  fit2 <- with(X, ash(betahat = beta_hat_2, sebetahat = seb2,
                      mixcompdist = "normal", prior="nullbiased",
                      optmethod = optmethod))

  sigma1 <- with(fit1$fitted_g, sd[!zapsmall(pi)==0])
  sigma2 <- with(fit2$fitted_g, sd[!zapsmall(pi)==0])
  mix_grid <- expand.grid("S1"=sigma1, "S2"=sigma2, "pi"=0)
  new_cause_grid(mix_grid)
}

new_cause_grid <- function(x = data.frame()){
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c("S1", "S2", "pi") %in% names(x)))
  structure(x, class = c("cause_grid", "data.frame"))
}


