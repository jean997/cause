#'@title Recompute elpd table for a CAUSE fit that is already computed
#'@param res A cause fit from running the cause function
#'@return A model table equivalent to res$elpd
#'@export
recompute_elpd_table <- function(res){
  mods <- res$elpd
  loos <- res$loos
  X <- new_cause_data(res$data)
  llmat_null <- loglik_loo(0, 0, 0,
                      res$sharing$rho, res$sharing$mix_grid$S1,
                      res$sharing$mix_grid$S2, res$sharing$mix_grid$pi,
                      X$beta_hat_1, X$beta_hat_2, X$seb1, X$seb2)


  for(j in seq(nrow(mods))){
    i1 <- which(names(loos) == mods$model1[j])
    i2 <- which(names(loos) == mods$model2[j])
    if(is.null(loos[[i1]])){
      diff <- llmat_null - loos[[i2]]$pointwise[,1]
      est <- sum(diff)
      se <- sd(diff)*sqrt(nrow(X))
    }else if(is.null(loos[i2])){
      diff <- loos[[i1]]$pointwise[,1] - llmat_null
      est <- sum(diff)
      se <- sd(diff)*sqrt(nrow(X))
    }else{
      comp <- loo_compare(loos[[i1]], loos[[i2]])
      if(rownames(comp)[1]=="model2"){
        est <- comp[2,1]
        se <- comp[2,2]
      }else{
        est <- -1*comp[2,1]
        se <- comp[2,2]
      }
    }
    mods$delta_elpd[j] <- est
    mods$se_delta_elpd[j] <- se
  }
  mods <- mods %>%
          mutate(z = delta_elpd/se_delta_elpd) %>%
          select(model1, model2, delta_elpd, se_delta_elpd, z)
  class(mods) <- c("cause_elpd", "data.frame")
  return(mods)
}

#'@title Estimate delta elpd for pairs of models
#'@param X Data
#'@param fits List of models of class cause_post.
#'@param variants Optional vector of SNPs to use for ELPD calculation
#'@param nsamps Number of samples to take from the posterior
#'@return A data frame with one row per pair of models giving the estimated
#'difference in elpds and standard error of that difference. If delta_elpd is negative
#'this is in favor of model 2 in the table. A positive delta_elpd is in favor of model 1.
#'@export
in_sample_elpd_loo <- function(X, fits, variants, nsamps=1000){

  stopifnot(inherits(X, "cause_data"))
  if(!missing(variants)){
    if(!all(variants %in% X$snp)){
      stop("Not all `variants` are in data.", call.=FALSE)
    }
    X <- filter(X, snp %in% variants)
    X <- new_cause_data(X)
  }
  stopifnot(length(fits) >=2)

  k <- length(fits)
  for(x in fits) stopifnot(inherits(x, "cause_post"))
  mods <- expand.grid(m1 = 1:k, m2=1:k)
  mods <- mods %>% filter(m1 < m2) %>%
                   mutate(delta_elpd = NA, se_delta_elpd = NA)
  if(!is.null(names(fits))){
    mods$model1 <- names(fits)[mods$m1]
    mods$model2 <- names(fits)[mods$m2]
  }else{
    mods$model1 <- m1
    mods$model2 <- m2
  }
  #For each data point, calculate log likelihood under posterior
  lls <- lapply(fits, FUN=function(fit){
               if(is.null(fit$joint_post)){
                 llmat <- loglik_loo(0, 0, 0,
                                     fit$rho, fit$mix_grid$S1,
                                     fit$mix_grid$S2, fit$mix_grid$pi,
                                     X$beta_hat_1, X$beta_hat_2, X$seb1, X$seb2)
                 return(llmat)
               }
               samps <-samp_from_grid(fit$joint_post, c("q", "gamma", "eta"), nsamps)
               llmat <- t(loglik_loo(samps$gamma, samps$eta, samps$q, fit$rho,
                                     fit$mix_grid$S1, fit$mix_grid$S2, fit$mix_grid$pi,
                                     X$beta_hat_1, X$beta_hat_2, X$seb1, X$seb2))
               return(llmat)
              })

  loos <- lapply(lls, function(llmat){
    if(ncol(llmat) == 1) return(NULL)
    loo(llmat, r_eff = rep(1, nrow(X)))
  })

  for(j in seq(nrow(mods))){
    i1 <- mods$m1[j]
    i2 <- mods$m2[j]
    if(is.null(loos[[i1]])){
      diff <- lls[[i1]] - loos[[i2]]$pointwise[,1]
      est <- sum(diff)
      se <- sd(diff)*sqrt(nrow(X))
    }else if(is.null(loos[i2])){
      diff <- loos[[i1]]$pointwise[,1] - lls[[i2]]
      est <- sum(diff)
      se <- sd(diff)*sqrt(nrow(X))
    }else{
      comp <- loo_compare(loos[[i1]], loos[[i2]])
      if(rownames(comp)[1]=="model2"){
        est <- comp[2,1]
        se <- comp[2,2]
      }else{
        est <- -1*comp[2,1]
        se <- comp[2,2]
      }
    }
    mods$delta_elpd[j] <- est
    mods$se_delta_elpd[j] <- se
  }
  mods <- mods %>% mutate(z = delta_elpd/se_delta_elpd) %>%
          select(model1, model2, delta_elpd, se_delta_elpd, z)
  class(mods) <- c("cause_elpd", "data.frame")
  ret <- list(mods = mods, loos = loos)
  return(ret)
}
