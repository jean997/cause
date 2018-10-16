#'@title Sample from grid approximation of posterior distribuitions
#'@description Sample from the joint posterior of parameters
#'@param joint_post Data frame giving the joint posterior
#'@param params Vector of names of parameters to sample
#'@param nsamps Number of samples to take
#'@return A data.frame that is nsamps by length(params) of samples
#'@export
samp_from_grid <- function(joint_post, params, nsamps){
  #will look for variable names x, xstart, xstop where x is a parameter
  #will also look for log_post
  stopifnot(all(params %in% names(joint_post)))
  stopifnot(all(paste0(params, "start") %in% names(joint_post)))
  stopifnot(all(paste0(params, "stop") %in% names(joint_post)))
  stopifnot("log_post" %in% names(joint_post))

  probs <- exp(joint_post$log_post)
  bucket <- sample(seq(nrow(joint_post)), size=nsamps, prob = probs, replace=T)

  samps <- sapply(seq_along(params), FUN=function(i){
    offset <- runif(nsamps, 0, 1)
    st <- paste0(params[i], "start")
    sp <- paste0(params[i], "stop")
    width <- joint_post[, sp]-joint_post[,st]
    joint_post[bucket, st] + offset*width[bucket]
  })
  samps <- data.frame(samps)
  names(samps) <- params
  return(samps)
}
