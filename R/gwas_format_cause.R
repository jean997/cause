#'@title Format data for CAUSE
#'@description Format GWAS summary statistics for CAUSE
#'@param X1 data.frame with data for GWAS 1
#'@param X2 data.frame with data for GWAS 2
#'@param snp_name_cols A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively. This is the column on which the data will be merged.
#'@param beta_hat_col A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively. If effect sizes are provided as odds ratios they must be converted
#'back into coeffecient estimates by taking the log.
#'@param se_cols A vector of length 2 specifying the name of the snp column
#'in X1 and X2 respectively.
#'@param A1_cols Column names for effect allele
#'@param A2_cols Column names for other allele
#'@details This function will try to merge data sets X1 and X2 on the specified columns. Where
#'necessary, it will flip the sign of effects so that the effect allele is the same in both
#'data sets. It will remove variants with ambiguous alleles or where the alleles (G/C or A/T) or
#'with alleles that do not match between data sets (e.g A/G in one data set and A/C in the other).
#'It will not remove variants that are simply strand flipped between the two data sets (e. g. A/C in one data set, T/G in the other).
#'@return An object of class cause_data and data.frame.
#'@export
gwas_format_cause <- function(X1, X2, snp_name_cols=c("snp", "snp"),
                              beta_hat_cols = c("beta_hat", "beta_hat"),
                              se_cols = c("se", "se"),
                              A1_cols= c("A1", "A1"),
                              A2_cols = c("A2", "A2")){

  X1 <- X1 %>% rename(snp = snp_name_cols[1],
                      beta_hat_1 = beta_hat_cols[1],
                      seb1 = se_cols[1],
                      A1 = A1_cols[1],
                      A2 = A2_cols[1]) %>%
        select(snp, beta_hat_1, seb1, A1, A2)

  X2 <- X2 %>% rename(snp = snp_name_cols[2],
                      beta_hat_2 = beta_hat_cols[2],
                      seb2 = se_cols[2],
                      A1 = A1_cols[2],
                      A2 = A2_cols[2]) %>%
        select(snp, beta_hat_2, seb2, A1, A2)
  cat("There are ", nrow(X1), " variants in GWAS 1 and ", nrow(X2),
      " variants in GWAS 2.\n")
  #Check case of allele coding
  upper1 <- X1$A1[1] == toupper(X1$A1[1])
  upper2 <- X2$A1[1] == toupper(X2$A1[1])

  if(upper1){
    if(!all(X1$A1 %in% c("A", "C", "T", "G"))){
      stop("Allele 1 in data set 1 contains illegal characters.")
    }
    if(!all(X1$A2 %in% c("A", "C", "T", "G"))){
      stop("Allele 2 in data set 1 contains illegal characters.")
    }
  }else{
    if(!all(X1$A1 %in% c("a", "c", "t", "g"))){
      stop("Allele 1 in data set 1 contains illegal characters.")
    }
    if(!all(X1$A2 %in% c("a", "c", "t", "g"))){
      stop("Allele 2 in data set 1 contains illegal characters.")
    }
  }
  if(upper2){
    if(!all(X2$A1 %in% c("A", "C", "T", "G"))){
      stop("Allele 1 in data set 2 contains illegal characters.")
    }
    if(!all(X2$A2 %in% c("A", "C", "T", "G"))){
      stop("Allele 2 in data set 2 contains illegal characters.")
    }
  }else{
    if(!all(X2$A1 %in% c("a", "c", "t", "g"))){
      stop("Allele 1 in data set 2 contains illegal characters.")
    }
    if(!all(X2$A2 %in% c("a", "c", "t", "g"))){
      stop("Allele 2 in data set 2 contains illegal characters.")
    }
  }

  X1 <- remove_ambiguous(X1, upper = upper1)
  X2 <- remove_ambiguous(X2, upper = upper2)

  X1 <- align_beta(X1, "beta_hat_1", upper1)
  X2 <- align_beta(X2, "beta_hat_2", upper2)

  cat("After removing ambiguous SNPs, there are ",
      nrow(X1), " variants in GWAS 1 and ", nrow(X2),
      " variants in GWAS 2.\n")


  X <-  inner_join(X1, X2, by="snp") %>%
        filter(!is.na(seb1) & !is.na(seb2) &
                 !is.na(beta_hat_1) & !is.na(beta_hat_2)) %>%
        filter(is.finite(seb1) & is.finite(seb2) &
                 is.finite(beta_hat_1) & is.finite(beta_hat_2)) %>%
        filter(seb1 > 0 & seb2 > 0)
  if(!upper1){
    X  <- X %>% mutate(A2.x = toupper(A2.x))
  }
  if(!upper2){
    X  <- X %>% mutate(A2.y = toupper(A2.y))
  }
  X <- filter(X, toupper(A2.x) == toupper(A2.y)) %>%
       select(-A1.y, -A2.y) %>%
       rename(A1 = A1.x, A2 = A2.x) %>%
      select(snp, beta_hat_1, seb1, beta_hat_2, seb2, A1, A2)
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
