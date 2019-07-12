#'@title Plot posteriors for one CAUSE fit
#'@param fit object of class cause_post `
#'@param type Either "posteriors" or "data". See details.
#'@param data If type=="data" then an object of type cause_data_fit must be supplied
#'@param pval_thresh p-value threshold used when type=="data".
#'@details If type == "posteriors", the function plots the marginal posterior
#'distribution of each parameter. If type=="data", the function plots data colored
#'by the probability that Z = 1.
#'@export
plot.cause_post <- function(fit, type=c("posteriors", "data"), data=NULL, pval_thresh=5e-8){
  type = match.arg(type)
  if(type == "posteriors"){
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
  if(type=="data"){
    stopifnot(inherits(data, "cause_data_fit"))
    if(all(c("gamma", "eta") %in% fit$params)){
      #Full model
      medians <- data.frame(param = c("gamma", "eta"), med = summary(fit)$quants[1, c(1, 2)])
      var <- "prob_Z1_full"
      title <- "Full Model"
    }else if("eta" %in% fit$params){
      #Confounding model
      medians <- data.frame(param = c( "eta"), med = summary(fit)$quants[1, 2])
      var <- "prob_Z1_conf"
      title <- "Confounding Only Model"
    }

    plt <- data %>% mutate(pval1 = 2*pnorm(-abs(beta_hat_1/seb1))) %>%
      filter(pval1 < pval_thresh) %>%
      ggplot(.) +
        geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
        geom_errorbar(aes(ymin = beta_hat_2 -1.96*seb2, ymax = beta_hat_2 + 1.96*seb2, x = beta_hat_1 ), color="grey") +
        geom_errorbarh(aes(y = beta_hat_2, xmin = beta_hat_1 - 1.96*seb1, xmax = beta_hat_1 + 1.96*seb1), color="grey") +
        geom_abline(aes(slope = med, linetype=param, intercept=0), data=medians) +
        geom_point(aes(x=beta_hat_1, y=beta_hat_2, col=get(var), size = -log10(pval1))) +
        scale_color_continuous(limits=c(0, 1), name = "P(Z = 1)", low="grey", high=muted("red")) +
        ggtitle(title) +
        theme_bw()
    return(plt)
    }
}

#'@title Plot CAUSE
#'@description Plot posteriors for confounding and causal models, display summary tables
#'@param res object of class cause
#'@param intern If TRUE, function returns a list of grobs. Otherwise it plots.
#'@param type Either "posteriors" or "data". See details.
#'@param pval_thresh p-value threshold used if type=="data".
#'@details If type == "posteriors", the function will plot the posterior distributions of the parameters
#'and display summary tables giving medians and credible intervals. If type == "data" the function
#'will plot the data thresholded on the trait 1 p-value if using pval_thresh.
#'@export
plot.cause <- function(res, intern=FALSE, type=c("posteriors", "data"), pval_thresh = 5e-8){
  type <- match.arg(type)
  if(type == "posteriors"){
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
    if(intern) return(plts)
    h <- arrangeGrob(grobs = plts,
                     layout_matrix = rbind(c(4, 4, 4, 3, 3, 3),
                                           c(NA, NA, 1, 1, 1, 1),
                                           c(2, 2, 2, 2,2,2)))
    plot(h)
  }
  if(type=="data"){
    plts <- lapply(c("conf", "full"), function(i){
      plt <- plot(res[[i]], type="data", data=res$data, pval_thresh=pval_thresh)
      return(plt)})
    max_delta_elpd <- max(abs(res$data$delta_elpd))
    plts[[3]] <-  res$data %>%
      mutate(pval1 = 2*pnorm(-abs(beta_hat_1/seb1))) %>%
      filter(pval1 < pval_thresh) %>%
      ggplot(.) +
        geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
        geom_errorbar(aes(ymin = beta_hat_2 -1.96*seb2, ymax = beta_hat_2 + 1.96*seb2, x = beta_hat_1 ), color="grey") +
        geom_errorbarh(aes(y = beta_hat_2, xmin = beta_hat_1 - 1.96*seb1, xmax = beta_hat_1 + 1.96*seb1), color="grey") +
        geom_point(aes(x=beta_hat_1, y=beta_hat_2, col=delta_elpd, size = -log10(pval1))) +
        scale_color_gradient2(name = "Contribution\nto test statisitc", mid = "grey", limits=c(-1, 1)*max_delta_elpd) +
        ggtitle("ELPD Contribution") +
        theme_bw()
    if(intern) return(plts)
    h <- arrangeGrob(grobs = plts, nrow=1)
    plot(h)
  }
}
