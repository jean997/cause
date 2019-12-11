
#'@title Estimate pi and mixture proportions under a specified model
#'@description Estimare the MAP for rho and the mixing proportions using
#'coordinate descent. Causal effect (gamma) and the effect of U (eta) are fixed at zero.
#'@param X An object of class cause_data containing data for the two traits.
#'@param mix_grid An object of class cause_grid containing variance pair candidates
#'@param rho_start Starting value for rho
#'@param z_prior_func Prior function for z = arctanh(rho)
#'@param null_wt Specifies the prior weight on the first entry of grid
#'@export
map_pi_rho_nonnull_fixed <- function(X, mix_grid, gamma, eta, q,
                                     rho_start=0, tol=1e-5, n.iter=20, null_wt = 10,
                                     z_prior_func = function(z){ dnorm(z, 0, 0.5, log=TRUE)},
			                               optmethod = c("mixSQP", "mixIP")){

  stopifnot(inherits(X, "cause_data"))
  stopifnot(inherits(mix_grid, "cause_grid"))
  optmethod <- match.arg(optmethod)
  if(optmethod=="mixSQP") optfun <- ashr:::mixSQP
    else optfun <- ashr:::mixIP
  if(!any(mix_grid$S1 == 0 & mix_grid$S2==0)){
    stop("Grid is invalid. It doesn't contain (0, 0).")
  }
  K <- nrow(mix_grid)
  p <- nrow(X)

  rho <- rho_old <- rho_start

  #Functions
  update_pi <- function(rho, gamma, eta, q){
    matrix_llik1 <- loglik_mat(rho, gamma, eta, q,
                              mix_grid$S1, mix_grid$S2,
                              X$beta_hat_1,
                              X$beta_hat_2,
                              X$seb1,
                              X$seb2)
    matrix_llik <- matrix_llik1 - apply(matrix_llik1, 1, max)
    matrix_lik <- exp(matrix_llik)
    w_res <- optfun(matrix_lik =matrix_lik,
                          prior=c(null_wt, rep(1, K-1)),
                          weights=rep(1, nrow(matrix_lik)))
    return(pmax(w_res$pihat, 0))
  }
  ddirichlet1 <- function(x, alpha, minx = 1e-6) {
    x <- pmax(x, minx)
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s <- sum((alpha - 1) * log(x))
    return(sum(s) - logD)
  }
  pi_prior_func <- function(pi){
	  ddirichlet1(pi, c(null_wt, rep(1, K-1)))
  	}
  arctanh <- function(rho){
    0.5*log((1+rho)/(1-rho))
  }
  li_func <- function(rho, pi){
    z <- arctanh(rho)
    ll <- loglik(rho, gamma, eta, q,
                    mix_grid$S1, mix_grid$S2,
                    pi,X$beta_hat_1,
                    X$beta_hat_2,
                    X$seb1,
                    X$seb2) +
      	   z_prior_func(z) +
	         pi_prior_func(pi)
     return(ll)

  }
  #rho
  li_func_rho <- function(rho){
    ll <- li_func(rho, pi)
    return(-ll)
  }

  #If there is no initial grid estimate
  cat("Initializing pi\n")
  if(all(mix_grid$pi==0)){
    pi <- pi_old <- update_pi(rho, gamma, eta, q)
  }else{
    pi <- pi_old <- mix_grid$pi
  }

  converged <- FALSE
  PIS <- matrix(pi, nrow=K, ncol=1)
  RHO <- c(rho)
  GAMMA <- c(gamma)
  ETA <- c(eta)
  Q <- c(q)
  ll <- li_func(rho, pi)
  LLS <- c(-1*ll)
  ct <-1

  while(!converged & ct <= n.iter){
    #Update rho
    cat("rho: ")
    opt_rho <-  optimize(f = li_func_rho, lower=-1, upper = 1, maximum=FALSE)
    rho <- opt_rho$minimum
    RHO <- c(RHO, rho)
    cat(rho, "\n")

    LLS <- c(LLS, -opt_rho$objective)

    #Update pi
    pi <- update_pi(rho, gamma, eta, q)
    ll <- li_func(rho, pi)

    LLS <- c(LLS, ll)
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
	      "gamma" = gamma, "eta" = eta, "q" = q,
              "loglik"=LLS[length(LLS)],
              "PIS"=PIS, "RHO"= RHO, "LLS"=LLS,
              "converged" = converged)
  class(fit) <- "cause_params"
  return(fit)

}

