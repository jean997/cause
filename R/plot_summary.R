
#'@title Summarize CAUSE Results
#'@description Summarize p-value testing that the causal model fits the data better and the posterior
#'medians and credible intervals for parameters
#'@param res object of class cause.
#'@param ci_size size of posterior credible intervals.
#'@param digits significant digits to report.
#'@export
summary.cause <- function(res, ci_size=0.95, digits=2){
  ans <- list()
  ans$quants <- lapply(c("conf", "full"), function(i){
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
  ans$z <- with(res$elpd, z[model1=="conf" & model2=="full"])
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
  ans$tab <- cbind(model=c("Confounding Only", "Causal"), ans$tab)
  ans$ci_size <- ci_size
  class(ans) <- "cause_summary"
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

#'@title Plot posteriors for one CAUSE fit
#'@param fit object of class cause_post `
#'@export
plot.cause_post <- function(fit){
  marge_post <- map_df(seq_along(fit$marge_post), function(i){
    p <- fit$marge_post[[i]]
    p$param <- fit$params[i]
    return(p)})
  marge_post <- select(marge_post, mid, width, param, post, prior) %>%
    gather("dist", "pdf", -mid, -param, -width)
  marge_post$param <- factor(marge_post$param, levels =c("gamma", "eta", "q"))
  plt <- ggplot(marge_post) + geom_line(aes(x=mid, y=pdf/width, linetype=dist)) +
          xlab("Parameter Value") + ylab("Density") +
          theme_bw() + theme(legend.title = element_blank()) +
          facet_wrap(~param, scale="free")
  return(plt)
}

#'@title Plot CAUSE
#'@description Plot posteriors for confounding and causal models, display summary tables
#'@param res object of class cause
#'@export
plot.cause <- function(res){
  plts <- lapply(c("conf", "full"), function(i){
    plt <- plot(res[[i]])
    return(plt)})
  elpd <- res$elpd %>% mutate(p = pnorm(z),
                              delta_elpd = signif(delta_elpd, digits=2),
                              se_delta_elpd = signif(se_delta_elpd, digits=2),
                              z = signif(z, digits=2),
                              p = signif(p, digits=2))
  elpd <- tableGrob(elpd)
  tab <- tableGrob(summary(res)$tab)
  plts[[3]] <- tab
  plts[[4]] <- elpd
  h <- arrangeGrob(grobs = plts,
                   layout_matrix = rbind(c(4, 4, 4, 3, 3, 3),
                                         c(NA, NA, 1, 1, 1, 1),
                                         c(2, 2, 2, 2,2,2)))
  plot(h)
}
