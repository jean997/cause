

#'@title Prior variance for gamma and eta
#'@description Choose prior variance for eta and gamma based on the data
#'@param X An object of class cause_data containing data for the two traits.
#'@param variants A vector of variants to include in the analysis.
#'@param prob See Details
#'@param pval_thresh A p-value threshold applied to trait 1 (M). This prevents
# inclusion of variants with beta_hat_1 very close to 0.
#'@details The function will return a variance, sigma, such that
#'P(abs(z) > k) = prob where z is a N(0, sigma) random variable, k = max(abs(beta_hat_2/beta_hat_1))
#'and prob is specified in the arguments of the function.
#'@return A prior variance for gamma and eta (scalar)
eta_gamma_prior <- function(X, variants, prob = 0.05, pval_thresh = 1e-3){
  stopifnot(inherits(X, "cause_data"))
  if("p1" %in% names(X)){
    X <- X %>% mutate(pval_m = p1)
  }else{
    X <- X %>% mutate(pval_m = 2*pnorm(-abs(beta_hat_1/seb1)))
  }
  if(!missing(variants)){
    X <- X %>% dplyr::filter(snp %in% variants)
  }
  X <- X %>%
    dplyr::filter(pval_m < pval_thresh) %>%
    new_cause_data()
  z <- with(X, max(abs(beta_hat_2/beta_hat_1), na.rm=TRUE))
  f <- function(sd){
    abs(prob/2 - pnorm(-z, sd=sd))
  }
  x <- optimize(f,lower = 0, upper = z)
  return(x$minimum)
}
