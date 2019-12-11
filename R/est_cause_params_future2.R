
#'@title Estimate CAUSE Nuisance Parameters
#'@description Estimates bivariate distribution of summary statistics and rho, the
#'correlation between summary statistics due to overlapping samples or population structure.
#'@param X An object of class cause_data containing data for the two traits.
#'@param variants A vector of variants to include. This list should be approximately LD
#'pruned and include variants genome wide.
#'@return An object of class cause_params
#'@export
est_cause_params_future2 <- function(X, variants,
				                            optmethod = c("mixSQP", "mixIP"), sigma_g_pval = 1e-3,
                                    null_wt = 10, max_candidates=Inf){
  optmethod <- match.arg(optmethod)
  stopifnot(inherits(X, "cause_data"))

  if(!all(variants %in% X$snp)){
    warning("Warning: Not all `variants` are in data.", call.=FALSE)
  }


  X <- filter(X, snp %in% variants)
  X <- new_cause_data(X)
  if(nrow(X) < 1e5){
    warning("Fewer than 100,000 variants are being used to estimate parametrs. ",
            "This can cause problems and is not recomended. You are using ", nrow(X),
            " variants.\n")
  }
  cat("Estimating CAUSE parameters with ", nrow(X), " variants.\n")
  cat("Running ashr to get variance pair candidates.\n")
  mix_grid <- variance_pair_candidates(X, optmethod=optmethod, max_candidates=max_candidates)
  cat("Done with that.\n")

  cat("Fitting the modal estimator\n")
  mr_dat <-try( X %>% filter(p_value < 5e-8 & ld_prune == TRUE) %>%
                         with(., MendelianRandomization::mr_input(bx = beta_hat_1, bxse = seb1,
                         by = beta_hat_2, byse = seb2,
                         snps = snp)), silent=TRUE)
  mbe_fit <- try(MendelianRandomization::mr_mbe(mr_dat, weighting="weighted", stderror="delta", phi=1), silent=TRUE)

  if(class(mbe_fit) == "try-error"){
    est <- 0
  }else{
    est <- mbe_fit@Estimate
  }

  params <- map_pi_rho_nonnull_fixed(X, mix_grid, gamma = est, eta = 0, q=0)

  #Filter out grid points with low mixing proportion
  params$mix_grid <- dplyr::filter(params$mix_grid, zapsmall(pi) > 0)

  return(params)
}


