
#'@title Posterior probability of acting through U
#'@description For a set of data, X, calculate the probability that each variant
#'acts through the shared factor U
#'@param X object of class cause_data. This need not be the same data used to compute 'fit'.
#'@param fit object of class cause_post.
#'@param nsamps Number of samples to take from the posterior (see Details).
#'@details
#'Let x_i represent the the data in the ith row of 'X'. Let Z_i be the indicator that
#'variant i affects U. We are interested in computing P(Z_i = 1 | posteriors) where
#'posteriors refers to the posterior parameter distributions contained in  'fit'.
#' We compute l(x_i | posteriors, Z=1) and l(x_i | posteriors, Z_i = 0) where l denotes the likelihood.
#' We then compute P(Z_i = 1 | posteriors) = l(x_i | posteriors, Z_i = 1)/(l(x_i | posteriors, Z_i = 1) + l(x_i | posteriors, Z_i = 0))
#' @return A vector of probabilities corresponding to the rows of X
Z_post <- function(X, fit, nsamps=1000, name_post = ""){
  if(!inherits(X, "cause_data")){
    warning("Converting X to cause_data.")
    X <- new_cause_data(X)
  }
  stopifnot(inherits(fit, "cause_post"))
  samps <-samp_from_grid(fit$joint_post,
                         c("q", "gamma", "eta"), nsamps)
  mix_grid <- fit$mix_grid
  rho <- fit$rho
  ll_mat_Z1 <- loglik_samps_Z1(samps$gamma, samps$eta, samps$q, rho,
                      mix_grid$S1, mix_grid$S2, mix_grid$pi,
                      X$beta_hat_1, X$beta_hat_2,
                      X$seb1, X$seb2)

  ll_mat_Z0 <- loglik_samps_Z0(samps$gamma, samps$eta, samps$q, rho,
                              mix_grid$S1, mix_grid$S2, mix_grid$pi,
                              X$beta_hat_1, X$beta_hat_2,
                              X$seb1, X$seb2)
  ll_Z1 = rowMeans(ll_mat_Z1)
  ll_Z0 = rowMeans(ll_mat_Z0)
  prob_Z1 <- exp(ll_Z1)/(exp(ll_Z1) + exp(ll_Z0))
  return(prob_Z1)
}
