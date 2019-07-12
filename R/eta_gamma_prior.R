

#'@title eta/gamma prior
#'@description Choose prior variance for eta and gamma empirically from data
#'@param X An object of class cause_data containing data for the two traits.
#'@param variants A vector of variants to include in the analysis.
#'@param prob See Details
#'@details The function will return a variance, sigma, such that
#'P(abs(z) > k) = prob where z is a N(0, sigma) random variable, k = max(abs(beta_hat_2/beta_hat_1))
#'and prob is specified in the arguments of the function.
#'@return A list with items conf, full, elpd, summary, and plot.
#'@export
eta_gamma_prior <- function(X, variants, prob = 0.05){
  stopifnot(inherits(X, "cause_data"))
  if(! missing(variants)){
    X <- filter(X, snp %in% variants)
  }
  z <- with(X, max(abs(beta_hat_2/beta_hat_1), na.rm=TRUE))
  f <- function(sd){
    abs(prob/2 - pnorm(-z, sd=sd))
  }
  x <- optimize(f,lower = 0, upper = z)
  return(x$minimum)
}
