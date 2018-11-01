
#'@title CAUSE posterior density estimation
#'@description Uses an adaptive grid approximation to estimate the posterior distribution of q, gamma, and eta or a subset of those
#'@param X An object of class cause_data
#'@param param_ests param_ests Object of class cause_params
#'@param params Vector of parameter names. Should be some subest of  c("q", "gamma", "eta") given in any order.
#'@param priors List of prior functions corresponding to parameters in params. Each element should be a
#'funciton that takes one argument and returns a log likelihood.
#'@param n_start Vector of length length(params) giving the number of starting bins for each dimension
#'@param max_post_per_bin Maximum posterior probability for each grid cube.
#'This parameter controls the fineness of the approximation.
#'@return An object of class cause_post.
#'@export
cause_grid_adapt <-  function(X, param_ests,
                              params = c("q", "gamma", "eta"),
                              priors = list(function(q){dbeta(q, 1, 10)},
                                            function(b){dnorm(b, 0, 0.6)},
                                            function(b){dnorm(b, 0, 0.6)}),
                              n_start = c(10, 20, 20),
                              max_post_per_bin = 0.001){

  stopifnot(inherits(X, "cause_data"))
  stopifnot(inherits(param_ests, "cause_params"))
  stopifnot(all(params %in% c("q", "gamma", "eta")))
  stopifnot(length(params)==length(priors))
  stopifnot(length(params)==length(n_start))

  ranges <- lapply(params, function(x){
    if(x == "q") return(c(0, 1))
    return(c(-1, 1))
  })
  range_fixed <- sapply(params, function(x){
    if(x=="q") return(TRUE)
    return(FALSE)
  })
  joint_post <- adapt2_grid(params, ranges, priors, n_start,
                            range_fixed,
                            param_ests$mix_grid, param_ests$rho,
                            X, max_post_per_bin)
  ranges <- joint_post$ranges
  joint_post <- joint_post$res
  post_marge <- marge_dists(joint_post, params, priors, ranges)

  R <- structure(list(joint_post = joint_post, marge_post = post_marge, params = params,
                      mix_grid = param_ests$mix_grid, rho = param_ests$rho, ranges=ranges, n_start = n_start,
                      max_post_per_bin = max_post_per_bin), class="cause_post")
  return(R)

}

#'@title Adaptive grid estimation of posteriors
#'@description Workhorse for posterior adaptive grid approximation. Called from cause_grid_adapt.
adapt2_grid <- function(params, ranges, priors, n_start,
                        range_fixed,
                        mix_grid, rho,
                        X, max_post_per_bin){
                        #end_bin_thresh=1e-8){

  k <- length(params)
  stopifnot(length(ranges) == k)
  stopifnot(length(priors) == k)
  stopifnot(length(n_start)==k)
  stopifnot(length(range_fixed)==k)

  vals <- get_vals2(ranges, n_start, priors)
  res <- post_from_vals(vals, params, X, mix_grid, rho)
  post_marge <- marge_dists(res, params, priors, ranges)
  range_set <- FALSE
  widths <- sapply(seq_along(params), function(i){(ranges[[i]][2]-ranges[[i]][1])/n_start[i]})
  n_add <- ceiling(n_start/2)
  end_bin_thresh <- sapply(post_marge, function(x){1/(nrow(x)*1000)})
  #Set range first then refine the grid
  while(!range_set){
    range_set <- TRUE
    for(i in seq_along(params)){
      if(range_fixed[i]) next
      n_new <- n_start
      n_new[i] <- n_add[i]
      if(with(post_marge[[i]], post[which.min(mid)]) > end_bin_thresh[i]){
        new_ranges <- ranges
        new_ranges[[i]][1] <- ranges[[i]][1]- n_add[i]*widths[i]
        new_ranges[[i]][2] <- ranges[[i]][1]
        ranges[[i]][1] <- new_ranges[[i]][1]

        vals <- get_vals2(new_ranges, n_new, priors)
        new_res <- post_from_vals(vals, params, X, mix_grid, rho)
        res <- rbind(res, new_res)
        res <- res %>% mutate(log_post  =  log_lik + log_prior,
                              log_post = log_post - logSumExp(log_post))

        range_set <- FALSE
	      n_start[i] <- n_start[i] + n_add[i]
      }
      if(with(post_marge[[i]], post[which.max(mid)]) > end_bin_thresh[i] ){
        new_ranges <- ranges
        new_ranges[[i]][2] <- ranges[[i]][2] +  n_add[i]*widths[i]
        new_ranges[[i]][1] <- ranges[[i]][2]
        ranges[[i]][2] <- new_ranges[[i]][2]
        vals <- get_vals2(new_ranges, n_new, priors)
        new_res <- post_from_vals(vals, params, X, mix_grid, rho)
        res <- rbind(res, new_res)
        res <- res %>% mutate(log_post  =  log_lik + log_prior,
                              log_post = log_post - logSumExp(log_post))

        range_set <- FALSE
	      n_start[i] <- n_start[i] + n_add[i]
      }
    }
    post_marge <-marge_dists(res, params, priors, ranges)
    end_bin_thresh <- sapply(post_marge, function(x){1/(nrow(x)*1000)})
  }
  #Refine grid
  thresh <- log(max_post_per_bin)
  n <- 3
  while(any(res$log_post > thresh)){
    ix <- which(res$log_post > thresh)
    nbrs <- unique(unlist(get_neighbors(ix, res, params)))
    ix <- unique(c(ix, nbrs))
    nr <- map_df(ix, function(i){
      new_ranges <- list()
      for(j in seq_along(params)){
        new_ranges[[j]] <- as.numeric(res[i, paste0(params[j], c("start", "stop"))])
      }
      vals <- get_vals2(new_ranges, rep(n, length(params)), priors)
      new_res <-post_from_vals(vals, params, X, mix_grid, rho)
      return(new_res)
    })
    res <- res[-ix,]
    res <- rbind(res, nr)
    res <- res %>% mutate(log_post  =  log_lik + log_prior,
                          log_post = log_post -logSumExp(log_post),
                          norm_log_post = log_post - log_volume)

  }
  R <- list(res = res, ranges = ranges)
  return(R)
}

#Takes in vals as produced by get_vals2 and outputs a data frame with all combos of starts/stops for each parameter
#head(res, 3)
#Var1 Var2 Var3          q     qstart      qstop     qwidth     qprior        gamma  gammastart  gammastop gammawidth  gammaprior
#5129    3    1    1 0.08333333 0.06666667 0.10000000 0.03333333 0.15293339  0.006666667  0.00000000 0.01333333 0.01333333 0.008864654
#5118    1    3    2 0.01666667 0.00000000 0.03333333 0.03333333 0.28752861 -0.006666667 -0.01333333 0.00000000 0.01333333 0.008864654
#5199    1    2    2 0.21666667 0.20000000 0.23333333 0.03333333 0.03721802  0.060000000  0.05333333 0.06666667 0.01333333 0.008820988
#eta   etastart     etastop   etawidth   etaprior   log_lik  log_prior  log_post log_volume norm_log_post
#5129 -0.1666667 -0.2000000 -0.13333333 0.06666667 0.04262911 -133.2121  -9.758654 -6.930810  -10.42674      3.495925
#5118 -0.1000000 -0.1333333 -0.06666667 0.06666667 0.04369367 -133.8830  -9.102668 -6.945751  -10.42674      3.480985
#5199 -0.1000000 -0.1333333 -0.06666667 0.06666667 0.04369367 -131.8405 -11.152136 -6.952689  -10.42674      3.474047
post_from_vals <- function(vals, params, X, mix_grid, rho){
  ix <- lapply(vals, function(x){seq(nrow(x))})
  res <- expand.grid(ix)
  for(i in seq_along(params)){
    res[[params[i]]] <- vals[[i]]$mid[ res[[paste0("Var", i)]]  ]
    res[[paste0(params[i], "start")]] <- vals[[i]]$begin[ res[[paste0("Var", i)]]  ]
    res[[paste0(params[i], "stop")]] <- vals[[i]]$end[ res[[paste0("Var", i)]]  ]
    res[[paste0(params[i], "width")]] <- vals[[i]]$width[ res[[paste0("Var", i)]]  ]
    res[[paste0(params[i], "prior")]] <- vals[[i]]$prior[ res[[paste0("Var", i)]]  ]
  }
  full_params <- c("gamma", "q", "eta")
  for(p in full_params){
    if(! p %in% params){
      res[[p]] <- 0
      res[[paste0(p, "start")]] <- 0
      res[[paste0(p, "stop")]] <- 0
    }
  }
  res[["log_lik"]] <- apply(res[, c("gamma","eta", "q")], 1, FUN = function(x){
    loglik(rho, x[1], x[2], x[3],
        mix_grid$S1, mix_grid$S2, mix_grid$pi,
        X$beta_hat_1, X$beta_hat_2,
        X$seb1, X$seb2)
  })
  res$log_prior <- Reduce("+", lapply(params, function(p){log(res[[paste0(p, "prior")]])}))
  res$log_volume <- Reduce("+", lapply(params, function(p){log(res[[paste0(p, "width")]])}))
  res <- res %>% mutate(log_post  =  log_lik + log_prior,
                        log_post = log_post - logSumExp(log_post))
  res$norm_log_post <- with(res, log_post-log_volume)
  return(res)
}

#Vals is a list of data frames with starts, stops, priors in each dimension
#length(vals) == length(params)
#lapply(vals, function(x){x[1:2,]})
#[[1]]
#begin end width  mid     prior
#1   0.0 0.1   0.1 0.05 0.6513216
#2   0.1 0.2   0.1 0.15 0.2413043
#
#[[2]]
#begin   end width   mid       prior
#1 -1.00 -0.96  0.04 -0.98 0.007008939
#2 -0.96 -0.92  0.04 -0.94 0.007797581
#
#[[3]]
#begin  end width  mid      prior
#1  -1.0 -0.8   0.2 -0.9 0.04342087
#2  -0.8 -0.6   0.2 -0.7 0.06744403
get_vals2 <- function(ranges, n, priors){
  k <- length(ranges)
  stopifnot(length(n)==k)
  stopifnot(length(priors)==k)
  vals <- list()
  for(i in seq_along(n)){
    s <- seq(ranges[[i]][1], ranges[[i]][2], length.out=n[i]+1)
    vals[[i]] <- data.frame("begin"=s[-(n[i]+1)], "end"=s[-1]) %>%
      mutate( width = end-begin, mid = begin + (width/2))
    vals[[i]]$prior <- apply(vals[[i]][,c("begin", "end")], 1, function(x){
      integrate(f=function(xx){priors[[i]](xx)}, lower=x[1], upper=x[2])$val})
  }
  return(vals)
}

get_neighbors <- function(ix, post, params){
  ivls <- map(params, function(p){
    Intervals(post[, paste0(p, c("start", "stop"))])
  })
  nbrs <- map(ivls, function(x){
    interval_overlap(x[ix,], x)
  }) %>% purrr::transpose() %>%
    map(., function(x){Reduce(dplyr::intersect, x)})
  return(nbrs)
}


neighbor_difference <- function(post, params, ix, thresh){
  ivls <- lapply(params, function(p){
    Intervals(post[, paste0(p, c("start", "stop"))])
  })
  ivl_overlap <- lapply(ivls, function(x){
    unlist(interval_overlap(x[ix,], x))
  })
  nbrs <- Reduce(dplyr::intersect, ivl_overlap)
  d <- post$norm_log_post[ix]-post$norm_log_post[nbrs]
  return(data.frame(ix=nbrs, dist=d))
  if(any(abs(d) > thresh)){
    return(c(ix, nbrs[abs(d) > thresh]))
  }else{
    return(NULL)
  }
}
