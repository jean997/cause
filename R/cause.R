
#'@title CAUSE
#'@description Fit CAUSE confounding only and full (causal) models,
#'calculate ELPD test statistic and estimate posteriors
#'@param X An object of class cause_data containing data for the two traits.
#'@param param_ests Object of class cause_params output by est_cause_params.
#'This contains estimates of the mixing proportions and an estimate of rho,
#'the correlation in test statistics that is due to overlapping samples or population structure.
#'@param variants A vector of variants to include in the analysis.
#'@param sigma_g Parameter specifying the prior distribution of gamma and eta.
#'gamma ~ N(0, sigma_g), eta ~ N(0, sigma_g).
#'@param qalpha,qbeta Parameters defining the prior distribution of q.
#'q ~ Beta(qalpha, qbeta)
#'@param max_q Largest value of q to be allowed. If max_q < 1 then the prior will be truncated.
#'@param force If true, do not give an error if parameter estimates did not converge.
#'@return A list with items conf, full, elpd, summary, and plot.
#'@export
cause <- function(X, param_ests, variants = X$snp,
                  sigma_g, qalpha = 1, qbeta=10,
                  max_q = 1, force=FALSE){
  stopifnot(inherits(X, "cause_data"))
  stopifnot(inherits(param_ests, "cause_params"))
  if(!param_ests$converged){
    if(!force){
      stop("The parameter estimates you are using did not converge.
           If you are sure you want to use them anyway, rerun with force=TRUE\n")
    }else{
      warning("Parameter estimates are not converged but we are going forward anyway.\n")
    }
  }
  if(!all(variants %in% X$snp)){
    stop("Not all `variants` are in data.", call.=FALSE)
  }
  X <- filter(X, snp %in% variants)
  X <- new_cause_data(X)
  cat("Estimating CAUSE posteriors using ", nrow(X), " variants.\n")
  stopifnot(inherits(param_ests, "cause_params"))

  m <- pbeta(max_q, qalpha, qbeta)
  qprior <- function(q){
    ret <- rep(0, length(q))
    ret[q > 0 & q < max_q] <- dbeta(q[q > 0 & q < max_q], qalpha, qbeta)/m
    return(ret)
  }

  if(missing(sigma_g)){
    sigma_g <- eta_gamma_prior(X)
  }

  fit2 <- cause_grid_adapt(X, param_ests,
                           max_post_per_bin = 0.001,
                           params = c("eta", "q"),
                           priors = list(function(b){dnorm(b, 0, sigma_g)},
                                         qprior),
                           n_start = c(20, 10))
  fit3 <- cause_grid_adapt(X, param_ests,
                               max_post_per_bin = 0.001,
                               params = c("gamma", "eta", "q"),
                               priors = list(function(b){dnorm(b, 0, sigma_g)},
                                             function(b){dnorm(b, 0, sigma_g)},
                                             qprior),
                               n_start = c(20, 20, 10))
  fit0 <- structure(list("joint_post"=NULL, rho = param_ests$rho, mix_grid=param_ests$mix_grid), class="cause_post")
  fits <- list("null"=fit0, "conf"=fit2, "full" = fit3)
  elpd <- in_sample_elpd_loo(X, fits)

  X$delta_elpd <- with(elpd, loos[[2]]$pointwise[,1] - loos[[3]]$pointwise[,1])

  X$prob_Z1_conf <- prob_confounding(X, fit2)
  X$prob_Z1_full <- prob_confounding(X, fit3)
  class(X) <- c("cause_data_fit", "cause_data", "data.frame")

  res <- list("conf"=fit2, "full" = fit3, elpd=elpd$mods, "loos" = elpd$loos, "data"=X,
              "sigma_g" = sigma_g, "qalpha" = qalpha, "qbeta" = qbeta)
  class(res) <- "cause"
  return(res)
}

#'@title Create a new cause_data_fit object
#'@param x a data.frame that includes columns snp, beta_hat_1, seb1, beta_hat_2, seb2, delta_elpd,
#'prob_Z1_conf, prob_Z1_full in any order.
#'x may also contain other columns
#'@return and object of class cause_data and data.frame.
#'@export
new_cause_data_fit <- function(x = data.frame()){
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2",
                  "delta_elpd", "prob_Z1_conf", "prob_Z1_full") %in% names(x)))
  structure(x, class = c("cause_data_fit", "cause_data", "data.frame"))
}
