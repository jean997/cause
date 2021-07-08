#'@title Format data for CAUSE
#'@description Format GWAS summary statistics for CAUSE
#'@param X1 data.frame with data for GWAS 1
#'@param X2 data.frame with data for GWAS 2
#'@param snp_name_cols A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively. This is the column on which the data will be merged.
#'@param beta_hat_cols A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively. If effect sizes are provided as odds ratios they must be converted
#'back into coeffecient estimates by taking the log.
#'@param se_cols A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively.
#'@param A1_cols Column names for effect allele
#'@param A2_cols Column names for other allele
#'@param pval_cols Column names for p-values. Can be ommitted or either element can be NA.
#'@param compute_pvals If p-values are missing in one or both studies, should the function
#'compute them. If true, p-values will be computed as `2*pnorm(-abs(beta_hat/se))`.
#'@details This function will try to merge data sets X1 and X2 on the specified columns. Where
#'necessary, it will flip the sign of effects so that the effect allele is the same in both
#'data sets. It will remove variants with ambiguous alleles or where the alleles (G/C or A/T) or
#'with alleles that do not match between data sets (e.g A/G in one data set and A/C in the other).
#'It will not remove variants that are simply strand flipped between the two data sets (e. g. A/C in one data set, T/G in the other).
#'@return An object of class cause_data and data.frame.
#'@export
gwas_merge <- function(X1, X2,
                       X1_formatted=FALSE, X2_formatted = FALSE,
                       snp_name_cols=c("snp", "snp"),
                       beta_hat_cols = c("beta_hat", "beta_hat"),
                       se_cols = c("se", "se"),
                       A1_cols= c("A1", "A1"),
                       A2_cols = c("A2", "A2"),
                       pval_cols  = c(NA,NA),
                       compute_pvals = TRUE){

  if(!X1_formatted){
    cat("Formatting X1\n")
    X1 <- gwas_format(X1, snp = snp_name_cols[1],
                      beta_hat = beta_hat_cols[1], se = se_cols[1],
                      A1 = A1_cols[1], A2 = A2_cols[1],
                      p_value = pval_cols[1],
                      compute_pval = compute_pvals)
  }else{
    X1 <- validate_cause_data_single(X1)
  }
  if(!X2_formatted){
    cat("Formatting X2\n")
    X2 <- gwas_format(X2, snp = snp_name_cols[2],
                      beta_hat = beta_hat_cols[2], se = se_cols[2],
                      A1 = A1_cols[2], A2 = A2_cols[2],
                      p_value = pval_cols[2],
                      compute_pval = compute_pvals)
  }else{
    X2 <- validate_cause_data_single(X2)
  }

  X <-  X1 %>%
        select(snp, beta_hat, se, A1, A2, p_value) %>%
        rename(beta_hat_1 = beta_hat,
               seb1 = se,
               p1 = p_value) %>%
        inner_join(., X2, by="snp") %>%
        select(snp, beta_hat_1, seb1, p1, beta_hat, se, p_value, A1.x, A2.x, A1.y, A2.y) %>%
        rename(beta_hat_2 = beta_hat,
               seb2 = se,
               p2 = p_value) %>%
        filter(!is.na(seb1) & !is.na(seb2) &
              !is.na(beta_hat_1) & !is.na(beta_hat_2)) %>%
        filter(is.finite(seb1) & is.finite(seb2) &
               is.finite(beta_hat_1) & is.finite(beta_hat_2)) %>%
        filter(seb1 > 0 & seb2 > 0) %>%
        filter(A2.x == A2.y) %>%
        select(-A1.y, -A2.y) %>%
        rename(A1 = A1.x, A2 = A2.x) %>%
        select(snp, beta_hat_1, seb1, p1, beta_hat_2, seb2, A1, A2, p2)
  cat("After merging and removing variants with inconsistent alleles, ",
      "there are ", nrow(X),
      " variants that are present in both studies and can be used with CAUSE.\n")
  new_cause_data(X)

}

remove_ambiguous <- function(X, upper=TRUE){
  if(upper){
    X <- X %>% dplyr::filter(!(A1 == "G" & A2 == "C") &
                             !(A1 == "C" & A2 == "G") &
                             !(A1 == "A" & A2 == "T") &
                             !(A1 == "T" & A2 == "A"))
    return(X)
  }
  X <- X %>% filter(!(A1 == "g" & A2 == "c") &
                             !(A1 == "c" & A2 == "g") &
                             !(A1 == "a" & A2 == "t") &
                             !(A1 == "t" & A2 == "a"))
  return(X)
}

#Flip signs and strabds so that allele 1 is allways A
align_beta <- function(X, beta_hat_name, upper=TRUE){
  flp = c("A" = "T", "G" = "C", "T" = "A",
          "C" = "G", "a"  = "t", "t" = "a",
          "c" = "g", "g" = "c")
  if(upper){
   X <- X %>% mutate( flip_strand = A1 == "T" | A2 == "T")
  }else{
   X <- X %>% mutate( flip_strand = A1 == "t" | A2 == "t")
  }
  X <- X %>% mutate(A1flp = case_when(flip_strand ~ flp[A1],
                                      TRUE ~ A1),
                    A2flp = case_when(flip_strand ~ flp[A2],
                                      TRUE ~ A2),
                    temp = case_when(A1flp == "A" | A1flp == "a" ~ get(beta_hat_name),
                                                   TRUE ~ -1*get(beta_hat_name))) %>%
            select(-A1, -A2) %>%
            mutate(A1 = case_when(A1flp == "A" | A1flp=="a" ~ A1flp,
                                   TRUE ~ A2flp),
                    A2 = case_when(A1flp == "A" | A1flp=="a" ~ A2flp,
                                   TRUE ~ A1flp)) %>%
            select(-A1flp, -A2flp, -flip_strand)
  ix <- which(names(X)==beta_hat_name)
  X <- X[,-ix]
  ix <- which(names(X) == "temp")
  names(X)[ix] <- beta_hat_name
  return(X)
}

#'@title Create a new cause_data object
#'@param x a data.frame that includes columns snp, beta_hat_1, seb1, beta_hat_2, seb2 in any order.
#'x may also contain other columns
#'@return and object of class cause_data and data.frame.
#'@export
new_cause_data <- function(x = data.frame()){
  stopifnot(inherits(x, "data.frame"))
  x <- validate_cause_data(x)
  #stopifnot(all(c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2") %in% names(x)))
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
    if(any(is.na(x[[n]]))){
      stop(paste0("Missing values in ", n), call.=FALSE)
    }
    if(any(!is.finite(x[[n]]))){
      stop(paste0("Infinite values in ", n), call.=FALSE)
    }
  }
  if(any(x$seb1 <= 0 | x$seb2 <=0)){
    stop("All of seb1 and seb2 must be positive", call.=FALSE)
  }
  x
}

validate_cause_data_single <- function(x){
  if(!inherits(x, "data.frame")){
    stop("`x` must inherit class data.frame", call.=FALSE)
  }
  req_names <- c("snp", "beta_hat", "se", "A1", "A2")
  if(!all(req_names %in% names(x))){
    stop(paste0("`x` must contain variables ", req_names),
         call.=FALSE)
  }
  for(n in c("beta_hat", "se")){
    if(any(is.na(x[[n]]))){
      stop(paste0("Missing values in ", n), call.=FALSE)
    }
    if(any(!is.finite(x[[n]]))){
      stop(paste0("Infinite values in ", n), call.=FALSE)
    }
  }
  if(any(x$se <= 0)){
    stop("All se must be positive", call.=FALSE)
  }
  if(any(x$A1 != "A")){
    stop("Alleles are not alligned. Use gwas_format to format data.", call.=FALSE)
  }
  x
}
