
#'@title Estimate pi and mixture proportions under null model
#'@description Estimare the MAP for rho and the mixing proportions using
#'coordinate descent. Causal and confounding effects are fixed at zero.
#'@param X An object of class cause_data containing data for the two traits.
#'@param mix_grid An object of class cause_grid containing variance pair candidates
#'@param rho_start Starting value for rho
#'@param z_prior_func Prior function for z = arctanh(rho)
#'@param null_wt Specifies the prior weight on the first entry of grid
#'@export
map_pi_rho <- function(X, mix_grid, rho_start=0,
                       tol=1e-7, n.iter=20, null_wt = 10,
                       z_prior_func = function(z){ dnorm(z, 0, 0.5, log=TRUE)}){

  stopifnot(inherits(X, "cause_data"))
  stopifnot(inherits(mix_grid, "cause_grid"))
  if(!any(mix_grid$S1 == 0 & mix_grid$S2==0)){
    stop("Grid is invalid. It doesn't contain (0, 0).")
  }
  K <- nrow(mix_grid)
  p <- nrow(X)

  rho <- rho_old <- rho_start
  #If there is no initial grid estimate
  if(all(mix_grid$pi==0)){
    matrix_llik1 <- loglik_mat(rho, 0, 0, 0,
                              mix_grid$S1, mix_grid$S2,
                              X$beta_hat_1,
                              X$beta_hat_2,
                              X$seb1,
                              X$seb2)
    matrix_llik = matrix_llik1 - apply(matrix_llik1, 1, max)
    matrix_lik = exp(matrix_llik)
    w_res = ashr:::mixIP(matrix_lik =matrix_lik,
                          prior=c(null_wt, rep(1, K-1)),
                          weights=rep(1, nrow(matrix_lik)))
    pi <- pi_old <- w_res$pihat
  }else{
    pi <- pi_old <- mix_grid$pi
  }
  pi_prior <- ddirichlet1(pi, c(null_wt, rep(1, K-1)))

  arctanh <- function(rho){
    0.5*log((1+rho)/(1-rho))
  }
  li_func <- function(rho){
    z <- arctanh(rho)
    loglik <- loglik(rho, 0, 0, 0,
                    mix_grid$S1, mix_grid$S2,
                    pi,X$beta_hat_1,
                    X$beta_hat_2,
                    X$seb1,
                    X$seb2) +
      z_prior_func(z) + pi_prior
    return(-loglik)
  }

  converged <- FALSE
  PIS <- matrix(pi, nrow=K, ncol=1)
  RHO <- c(rho)
  LLS <- c(-1*li_func(rho))
  ct <-1


  while(!converged & ct <= n.iter){
    #Update rho
    opt_rho <-  optimize(f = li_func, lower=-1, upper = 1, maximum=FALSE)
    rho <- opt_rho$minimum
    LLS <- c(LLS, -opt_rho$objective)
    RHO <- c(RHO, rho)
    #Update pi
    matrix_llik1 <- loglik_mat(rho, 0, 0, 0,
                              mix_grid$S1, mix_grid$S2,
                              X$beta_hat_1,
                              X$beta_hat_2,
                              X$seb1,
                              X$seb2)
    matrix_llik = matrix_llik1 - apply(matrix_llik1, 1, max)
    matrix_lik = exp(matrix_llik)
    w_res = ashr:::mixIP(matrix_lik =matrix_lik,
                         prior=c(null_wt, rep(1, K-1)),
                         weights=rep(1, nrow(matrix_lik)))
    pi <- w_res$pihat
    pi_prior <- ddirichlet1(pi, c(null_wt, rep(1, K-1)))
    loglik <- li_func(rho)

    LLS <- c(LLS, -1*loglik)
    PIS <- cbind(PIS, pi)

    #Test for convergence
    test <- max(abs(c(rho, pi)-c(rho_old, pi_old)))
    cat(ct, test, "\n")
    if(test < tol) converged <- TRUE
    rho_old <- rho
    pi_old <- pi
    ct <- ct + 1
  }
  if(!all(diff(LLS) > -1e-4)) cat("Warning: This may not be a local maximum ", min(diff(LLS)), "\n")
  mix_grid$pi <- pi
  fit <- list("rho"=rho, "pi"=pi, "mix_grid"=mix_grid,
              "loglik"=LLS[length(LLS)],
              "PIS"=PIS, "RHO"= RHO, "LLS"=LLS,
              "converged" = converged)
  fit$prior <- z_prior_func(arctanh(rho)) + pi_prior
  hes <- hessian(li_func, rho)
  fit$var <- solve(hes)
  class(fit) <- "cause_params"
  return(fit)

}

ddirichlet1 <- function(x, alpha) {
  logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  s <- sum((alpha - 1) * log(x))
  return(sum(s) - logD)
}
