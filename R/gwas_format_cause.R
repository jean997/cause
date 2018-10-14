#'@title Format data for CAUSE
#'@description Format GWAS summary statistics for CAUSE
#'@param X1 data.frame with data for GWAS 1
#'@param X2 data.frame with data for GWAS 2
#'@param snp_name_cols A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively. This is the column on which the data will be merged.
#'#'@param beta_hat_col A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively. If effect sizes are provided as odds ratios they must be converted
#'back into coeffecient estimates by taking the log.
#'#'@param se_cols A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively.
#'@return An object of class cause_data
#'@export
gwas_format_cause <- function(X1, X2, snp_name_cols=c("snp", "snp"),
                              beta_hat_cols = c("beta_hat", "beta_hat"),
                              se_cols = c("se", "se")){

  X1 <- X1 %>% rename(snp = snp_name_cols[1],
                      beta_hat_1 = beta_hat_cols[1],
                      seb1 = se_cols[1]) %>%
        select(snp, beta_hat_1, seb1)

  X2 <- X2 %>% rename(snp = snp_name_cols[2],
                      beta_hat_2 = beta_hat_cols[2],
                      seb2 = se_cols[2]) %>%
        select(snp, beta_hat_2, seb2)
  cat("There are ", nrow(X1), " variants in GWAS 1 and ", nrow(X2),
      " variants in GWAS2.\n")
  X <-  inner_join(X1, X2, by="snp") %>%
        filter(!is.na(seb1) & !is.na(seb2) &
                 !is.na(beta_hat_1) & !is.na(beta_hat_2)) %>%
        filter(is.finite(seb1) & is.finite(seb2) &
                 is.finite(beta_hat_1) & is.finite(beta_hat_2)) %>%
        filter(seb1 >= 0 & seb2 >= 0)

  cat("There are ", nrow(X), " variants overlapping which can be used with CAUSE.")
  new_cause_data(X)

}

new_cause_data <- function(x = data.frame()){
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2") %in% names(x)))
  structure(x, class = c("cause_data", "data.frame"))
}

validate_cause_data <- function(x){
  if(!inherits(x, "data.frame")){
    stop("`x` must inherit class data.frame", call.=FALSE)
  }
  req_names <- c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2")
  if(!all(req_names %in% names(x))){
    stop(paste0("`x` must contain variables ", req_names),
         call.=FALSE)
  }
  for(n in req_names[-1]){
    if(any(is.na(X[[n]]))){
      stop(paste0("Missing values in ", n), call.=FALSE)
    }
    if(any(!is.finite(X[[n]]))){
      stop(paste0("Infinite values in ", n), call.=FALSE)
    }
  }
  if(any(x$seb1 <= 0 | x$seb2 <=0)){
    stop("All of seb1 and seb2 must be positive", call.=FALSE)
  }
  x
}
