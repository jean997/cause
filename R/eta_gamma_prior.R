#'@export
eta_gamma_prior <- function(X, prob = 0.05){
  z <- with(X, max(abs(beta_hat_2/beta_hat_1), na.rm=TRUE))
  f <- function(sd){
    abs(prob/2 - pnorm(-z, sd=sd))
  }
  x <- optimize(f,lower = 0, upper = z)
  return(x$minimum)
}
