
#'@title Estimate pi and mixture proportions under a specified model
#'@description Estimare the MAP for rho and the mixing proportions using
#'coordinate descent. Causal effect (gamma) and the effect of U (eta) are fixed at zero.
#'@param X An object of class cause_data containing data for the two traits.
#'@param mix_grid An object of class cause_grid containing variance pair candidates
#'@param rho_start Starting value for rho
#'@param z_prior_func Prior function for z = arctanh(rho)
#'@param null_wt Specifies the prior weight on the first entry of grid
#'@export
map_pi_rho_nonnull_qge <- function(X, mix_grid, sigma_g, qalpha, qbeta,
			       q_start = qbeta(0.5, qalpha, qbeta),
			       gamma_start = 0, eta_start = 0,
             rho_start=0, tol=1e-5, n.iter=20, null_wt = 10,
             z_prior_func = function(z){ dnorm(z, 0, 0.5, log=TRUE)},
			       q_prior_func = function(q){dbeta(q, qalpha, qbeta, log=TRUE)},
			       gammaeta_prior_func = function(g){dnorm(g, 0, sigma_g, log=TRUE)},
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
  gamma <- gamma_old <- gamma_start
  eta <- eta_old <- eta_start
  q <- q_old <- q_start

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

  pi_prior_func <- function(pi){
	  ddirichlet1(pi, c(null_wt, rep(1, K-1)))
  	}
  arctanh <- function(rho){
    0.5*log((1+rho)/(1-rho))
  }
  li_func <- function(rho, gamma, eta, q, pi){
    z <- arctanh(rho)
    ll <- loglik(rho, gamma, eta, q,
                    mix_grid$S1, mix_grid$S2,
                    pi,X$beta_hat_1,
                    X$beta_hat_2,
                    X$seb1,
                    X$seb2) +
      	   z_prior_func(z) +
	   pi_prior_func(pi) +
	   gammaeta_prior_func(gamma)+
	   gammaeta_prior_func(eta) +
	   q_prior_func(q)
     return(ll)

  }
  #rho
  li_func_rho <- function(rho){
    ll <- li_func(rho, gamma, eta, q, pi)
    return(-ll)
  }
  #gamma
  li_func_gamma <- function(gamma){
    ll <- li_func(rho, gamma, eta, q, pi)
    return(-ll)
  }
  #eta
  li_func_eta <- function(eta){
    ll <- li_func(rho, gamma, eta, q, pi)
    return(-ll)
  }

  #q
  li_func_q <- function(q){
    ll <- li_func(rho, gamma, eta, q, pi)
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
  ll <- li_func(rho, gamma, eta, q, pi)
  LLS <- c(-1*ll)
  ct <-1

  while(!converged & ct <= n.iter){
    #Update rho
    cat("rho: ")
    opt_rho <-  optimize(f = li_func_rho, lower=-1, upper = 1, maximum=FALSE)
    rho <- opt_rho$minimum
    RHO <- c(RHO, rho)
    cat(rho)
    if(ct <= 2){
      #update q
      opt_q <-  optimize(f = li_func_q, lower=0, upper =1, maximum=FALSE)
      q <- opt_q$minimum
      Q <- c(Q, q)
      cat(" q: ", q)
      #update gamma
      opt_gamma <-  optimize(f = li_func_gamma, lower=-3*sigma_g, upper = 3*sigma_g, maximum=FALSE)
      gamma <- opt_gamma$minimum
      GAMMA <- c(GAMMA, gamma)
      cat(" gamma: ", gamma)
      #update eta
      opt_eta <-  optimize(f = li_func_eta, lower=-3*sigma_g, upper = 3*sigma_g, maximum=FALSE)
      eta <- opt_eta$minimum
      ETA <- c(ETA, eta)
      cat( " eta: ", eta, "\n")
    }
    ll <- li_func(rho, gamma, eta, q, pi)
    LLS <- c(LLS, ll)

    #Update pi
    pi <- update_pi(rho, gamma, eta, q)
    ll <- li_func(rho, gamma, eta, q, pi)

    LLS <- c(LLS, ll)
    PIS <- cbind(PIS, pi)

    #Test for convergence
    test <- max(abs(c(gamma, eta, q, rho, pi)-c(gamma_old, eta_old, q_old, rho_old, pi_old)))
    cat(ct, test, "\n")
    if(test < tol) converged <- TRUE
    gamma_old <- gamma
    eta_old <- eta
    q_old <- q
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

