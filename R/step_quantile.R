
#'@title Quantile of a step pdf
#'@param x quantile or vector of quantiles
#'@param start vector of step starts
#'@param stop vector of step stops
#'@param p vector of step probabilities
#'@export
step_quantile <- function(x, start, stop, p){
  stopifnot(length(start)==length(stop))
  stopifnot(length(stop)==length(p))
  cdf <- cumsum(p)
  q <- sapply(x, FUN=function(xx){
    ix <- min(which(cdf > xx))
    if(ix==1) res <- xx
      else res <- xx-cdf[ix-1]
    (res/p[ix])*(stop[ix]-start[ix]) + start[ix]
  })
  return(q)
}
