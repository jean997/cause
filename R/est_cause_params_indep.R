
#'@title Estimate CAUSE Nuisance Parameters using independence grid
#'@description Estimates bivariate distribution of summary statistics and rho, the
#'correlation between summary statistics due to overlapping samples or population structure.
#'@param X An object of class cause_data containing data for the two traits.
#'@param variants A vector of variants to include. This list should be approximately LD
#'pruned and include variants genome wide.
#'@return An object of class cause_params
#'@export
est_cause_params_indep <- function(X, variants, gamma, eta, q,
                                   optmethod = c("mixSQP", "mixIP"),
                                   z_prior_func = function(z){ dnorm(z, 0, 0.5, log=TRUE)},
                             null_wt = 10, gridmult=sqrt(2)){
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


  fit1 <- with(X, ash(betahat = beta_hat_1, sebetahat = seb1,
                      mixcompdist = "normal", prior="nullbiased",
                      optmethod=optmethod, gridmult=gridmult))


  fit2 <- with(X, ash(betahat = beta_hat_2, sebetahat = seb2,
                      mixcompdist = "normal", prior="nullbiased",
                      optmethod = optmethod, gridmult=gridmult))

  mix_grid <- expand.grid(S1 = seq_along(fit1$fitted_g$sd), S2 = seq_along(fit2$fitted_g$sd))
  mix_grid$pi <- fit1$fitted_g$pi[mix_grid$S1]*fit2$fitted_g$pi[mix_grid$S2]
  mix_grid <- filter(mix_grid, zapsmall(mix_grid$pi) > 0)
  mix_grid <- mutate(mix_grid, S1 = fit1$fitted_g$sd[S1],
                        S2 = fit1$fitted_g$sd[S2])


  arctanh <- function(rho){
    0.5*log((1+rho)/(1-rho))
  }
  li_func <- function(rho){
    z <- arctanh(rho)
    ll <- loglik(rho, gamma, eta, q,
                 mix_grid$S1, mix_grid$S2,
                 mix_grid$pi,X$beta_hat_1,
                 X$beta_hat_2,
                 X$seb1,
                 X$seb2) +
      z_prior_func(z)
    return(ll)

  }
  opt_rho <-  optimize(f = li_func, lower=-1, upper = 1, maximum=TRUE)
  rho <- opt_rho$maximum

  params <- list("rho" = rho, "pi"=mix_grid$pi, "mix_grid"=mix_grid, converged=TRUE)
  class(params) <- "cause_params"
  #params <- map_pi_rho(X, mix_grid, optmethod=optmethod, null_wt = null_wt)
  #Filter out grid points with low mixing proportion
  #params$mix_grid <- dplyr::filter(params$mix_grid, zapsmall(pi) > 0)

  return(params)
}



