#'@title Sample from grid approximation of posterior distribuitions
#'@description Sample from the joint posterior of parameters
#'@param joint_post Data frame giving the joint posterior
#'@param params Vector of names of parameters to sample
#'@params nsamps Number of samples to take
#'@retunr A data.frame that is nsamps by length(params) of samples
#'@export
samp_from_grid <- function(post, params, nsamps){
  #will look for variable names x, xstart, xstop where x is a parameter
  #will also look for log_post
  stopifnot(all(params %in% names(post)))
  stopifnot(all(paste0(params, "start") %in% names(post)))
  stopifnot(all(paste0(params, "stop") %in% names(post)))
  stopifnot("log_post" %in% names(post))

  probs <- exp(post$log_post)
  bucket <- sample(seq(nrow(post)), size=nsamps, prob = probs, replace=T)

  samps <- sapply(seq_along(params), FUN=function(i){
    offset <- runif(nsamps, 0, 1)
    st <- paste0(params[i], "start")
    sp <- paste0(params[i], "stop")
    width <- post[, sp]-post[,st]
    post[bucket, st] + offset*width[bucket]
  })
  samps <- data.frame(samps)
  names(samps) <- params
  return(samps)
}
