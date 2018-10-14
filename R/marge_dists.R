marge_dists <- function(joint_post, params, priors, ranges){
  post_marge <- list()
  for(i in seq_along(params)){
    param <- params[i]
    r1 <- ranges[[i]][1]
    r2 <- ranges[[i]][2]

    start_nm <- paste0(param, "start")
    stop_nm <- paste0(param, "stop")
    joint_post[[start_nm]] <- round(joint_post[[start_nm]], digits=6)
    joint_post[[stop_nm]] <- round(joint_post[[stop_nm]], digits=6)
    width_nm <- paste0(param, "width")

    min_width <- min(joint_post[,width_nm])

    starts <- seq(r1, r2, by=min_width) %>% round(., digits=6)
    stopifnot(all(unique(joint_post[[start_nm]]) %in% starts))
    post_marge[[i]] <- data.frame("begin"=starts[-length(starts)], "end"=starts[-1]) %>%
      mutate( width = end-begin, mid = begin + (width/2))
    post_marge[[i]]$prior <- apply(post_marge[[i]][,c("begin", "end")], 1, function(x){
      integrate(f=function(xx){priors[[i]](xx)}, lower=x[1], upper=x[2])$val})

    post_marge[[i]]$post <- apply(post_marge[[i]], 1, FUN=function(x){
      strt <- x[1]
      ix <- which( joint_post[,start_nm] <= strt & joint_post[,stop_nm] > strt)
      with(joint_post[ix,], exp(logSumExp(log_post - log(get(width_nm)) + log(min_width))))
    })
  }
  return(post_marge)
}
