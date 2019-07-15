
#'@title Summarize CAUSE Results
#'@description Summarize p-value testing that the causal model fits the data better and the posterior
#'medians and credible intervals for parameters
#'@param res object of class cause.
#'@param ci_size size of posterior credible intervals.
#'@param digits significant digits to report.
#'@export
summary.cause <- function(res, ci_size=0.95, digits=2){
  ans <- list()
  ans$quants <- lapply(c("sharing", "causal"), function(i){
    fit <- res[[i]]
    gamma <- eta <- q <- NA
    full_params <- c("gamma", "eta", "q")
    tt <- sapply(full_params, function(p){
      if(!p %in% fit$params){
        return(rep(NA, 3))
      }else{
        ix <- which(fit$params == p)
        qs <- with(fit$marge_post[[ix]], step_quantile(c(0.5, (1-ci_size)/2, 1-((1-ci_size)/2)),
                                                       begin, end, post))
        return(qs)
      }
    })
    return(tt)
  })
  ans$z <- with(res$elpd, z[model1=="sharing" & model2=="causal"])
  ans$p <- pnorm(ans$z)
  ans$tab <- sapply(ans$quants, function(y){
    apply(y, 2, function(z){
      if(is.na(z[1])) return(NA)
      paste0(round(z[1], digits=digits), " (",
             round(z[2], digits=digits), ", ",
             round(z[3], digits=digits), ")")
    })
  })
  ans$tab <- t(ans$tab)
  ans$tab <- cbind(model=c("Sharing", "Causal"), ans$tab)
  ans$ci_size <- ci_size
  class(ans) <- "cause_summary"
  return(ans)
}

#'@title Summarize CAUSE Results for a single fit
#'@description Summarize posterior
#'medians and credible intervals for parameters
#'@param fit object of class cause_post.
#'@param ci_size size of posterior credible intervals.
#'@param digits significant digits to report.
#'@export
summary.cause_post <- function(fit, ci_size = 0.95, digits=2){
  ans <- list()

  full_params <- c("gamma", "eta", "q")
  ans$quants <-  sapply(full_params, function(p){
      if(!p %in% fit$params){
        return(rep(NA, 3))
      }else{
        ix <- which(fit$params == p)
        qs <- with(fit$marge_post[[ix]], step_quantile(c(0.5, (1-ci_size)/2, 1-((1-ci_size)/2)),
                                                       begin, end, post))
        return(qs)
      }
    })
  ans$tab <- apply(ans$quants, 2, function(z){
      if(is.na(z[1])) return(NA)
      paste0(round(z[1], digits=digits), " (",
             round(z[2], digits=digits), ", ",
             round(z[3], digits=digits), ")")
    })
  ans$ci_size <- ci_size
  class(ans) <- "summary_cause_post"
  return(ans)

}


#'@title Print CAUSE Summary
#'@description Print a CAUSE summary
#'@param x object of class cause_summary
#'@param digits significant digits to report.
#'@export
print.cause_summary <- function(x, digits=2){
  cat("p-value testing that causal model is a better fit: ", signif(x$p, digits=digits), "\n")
  cat("Posterior medians and ", 100*x$ci_size, "% credible intervals:\n")
  print(x$tab, row.names=FALSE)
}

#'@title Print CAUSE fit Summary
#'@description Print a CAUSE fit summary
#'@param x object of class cause_summary_post
#'@export
print.summary_cause_post <- function(x){
  ix <- is.na(x$quants[1,])
  if(all(ix == FALSE)) cat("Full Model:\n")
    else if(ix[1] == TRUE & ix[2]==FALSE & ix[3] == FALSE) cat("Sharing Model:\n")
      else cat("Unrecognized Model:\n")
  cat("Posterior medians and ", 100*x$ci_size, "% credible intervals:\n")
  print(x$tab)
}
