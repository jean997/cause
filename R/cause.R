
#'@title CAUSE
#'@description Fit CAUSE confounding only and full (causal) models,
#'calculate ELPD test statistic and estimate posteriors
#'@param X An object of class cause_data containing data for the two traits.
#'@param param_ests Object of class cause_params output by est_cause_params.
#'This contains estimates of the mixing proportions and an estimate of rho,
#'the correlation in test statistics that is due to overlapping samples or population structure.
#'@param sigma_g Parameter specifying the prior distribution of gamma and eta.
#'gamma ~ N(0, sigma_g), eta ~ N(0, sigma_g).
#'@param qalpha,qbeta Parameters defining the prior distribution of q.
#'q ~ Beta(qalpha, qbeta)
#'@return A list with items conf, full, elpd, summary, and plot.
#'@export
cause <- function(X, variants, param_ests,
                  sigma_g = 0.6, qalpha = 1, qbeta=10){
  stopifnot(inherits(X, "cause_data"))
  if(!all(variants %in% X$snp)){
    stop("Not all `variants` are in data.", call.=FALSE)
  }
  X <- filter(X, snp %in% variants)
  X <- new_cause_data(X)
  cat("Estimating CAUSE posteriors using ", nrow(X), " variants.\n")
  stopifnot(inherits(param_ests, "cause_params"))

  fit2 <- cause_grid_adapt(X, param_ests,
                           max_post_per_bin = 0.001,
                           params = c("eta", "q"),
                           priors = list(function(b){dnorm(b, 0, sigma_g)},
                                             function(q){dbeta(q, qalpha, qbeta)}),
                           n_start = c(20, 10))
  fit3 <- cause_grid_adapt(X, param_ests,
                               max_post_per_bin = 0.001,
                               params = c("gamma", "eta", "q"),
                               priors = list(function(b){dnorm(b, 0, sigma_g)},
                                             function(b){dnorm(b, 0, sigma_g)},
                                             function(q){dbeta(q, qalpha, qbeta)}),
                               n_start = c(20, 20, 10))
  fit0 <- structure(list("joint_post"=NULL, rho = param_ests$rho, mix_grid=param_ests$mix_grid), class="cause_post")
  fits <- list("null"=fit0, "conf"=fit2, "full" = fit3)
  elpd <- in_sample_elpd_loo(X, fits)
  res <- list("conf"=fit2, "full" = fit3, elpd=elpd)
  class(res) <- "cause"
  return(res)
}
