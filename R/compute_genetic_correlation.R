
compute_gc_one <- function(gamma, eta, q, mix_grid){
  var_bm <- with(mix_grid, sum(pi*S1^2))
  var_theta <- with(mix_grid, sum(pi*S2^2))
  rho_g <- (gamma + q*eta)*var_bm/sqrt((gamma + q*eta)*(var_bm^2) + var_bm*var_theta)
  return(rho_g)
}

compute_genetic_correlation <- function(fit, nsamps=1000){
  samps <-samp_from_grid(fit$joint_post, c("q", "gamma", "eta"), nsamps)
  rhos <- apply(samps,1, function(x){compute_gc_one(x["gamma"], x["eta"], x["q"], fit$mix_grid)})
  return(rhos)
}

