
#'@title Prune variants for LD
#'@description Get an LD pruned list of variants optionally prioritizing variants with low p-values
#'@param variants Either a vector of SNP names or a data.frame. If a data.frame, provide pval_cols
#'containing the names of pvalue columns to use for pruning. Results will contain one list per p-value
#'column. If variants is a vector, LD pruning will be random.
#'@param ld data.frame with three columns: colsnp, rowsnp, and r2. Each row corresponds to an entry
#' of the LD matrix. If a pair of variants does not appear, they are assumed to have LD of 0 (or below the threshold).
#' @param total_ld_variants vector of variant names used to construct 'ld'.
#' @param pval_cols Vector of names of pval_cols. If not missing, variangs must be a data.frame.
#' pval_cols may contain NA entries. NA's indicate that you desire a pruned set that is not
#' prioritized by LD.
#' @param pval_thresh A vector the same length as 'pval_cols'. Only variants with p-value
#' below the threshold will be retained. Not required if variants is a vector. Entries of pval_thresh
#' may be Inf.
#' @param r2_thresh r^2 threshold for pruning. The final list will contain no pairs
#' with r2 > r2_thresh
#' @param variant_name If variants is a data.frame, the column name containing the variant id.
#' @details ld_prune will generate any desired number of LD pruned variant lists. At minimum, it requires
#' three arguments: 'variants' is the set of variants to be pruned, 'ld' is a data frame giving estimated
#' LD from a reference panel and 'total_ld_variants' is the list of variants used to construct the ld data.
#' Only variants ovlapping between 'variants' and 'total_ld_variants' will be considered.
#'
#' If only these three arguments are given, the returned object will be a vector of LD pruned variants.
#' These will be pruned randomly as no p-value information has been supplied.
#'
#' Another option is to prune variants prefentially retaining those with small p-values.
#' For this option, the 'variants' argument should be a data.frame that includes at minimum
#' a column containing the variant name and a column containing the pvalue. You can prune on
#' multiple p-values at once which can save time.
#'@export
ld_prune <- function(variants, ld, total_ld_variants, pval_cols, pval_thresh,
                     r2_thresh = 0.1,  variant_name="snp"){
  #check format of LD data frame
  if(!inherits(ld, "data.frame")){
    stop("ld should be a data frame with columns rowsnp, colsnp, r2")
  }
  if(!all(c("rowsnp", "colsnp", "r2") %in% names(ld))){
    stop("ld should be a data frame with columns rowsnp, colsnp, r2")
  }
  if(!missing(pval_cols)){
    stopifnot(inherits(variants, "data.frame"))
    if(missing(pval_thresh)) pval_thresh <- rep(Inf, length(pval_cols))
    else stopifnot(length(pval_cols)==length(pval_thresh))
  }
  if(is.vector(variants)){
    cat("Pruning for LD randomly (no p-values supplied)\n")
    n <- length(variants)
    variants <- data.frame(snp = variants, 
			   pval = sample(seq(n), size=n, replace=FALSE)/n,
			   stringsAsFactors=FALSE)
    pval_cols <- c("pval")
    variant_name <- "snp"
    pval_thresh <- Inf
  }

  stopifnot(inherits(variants, "data.frame"))
  stopifnot(variant_name %in% names(variants))
  n <- nrow(variants)
  cat("You have suppplied information for ", n, " variants.\n")

  if(any(is.na(pval_cols))){
    k <- sum(is.na(pval_cols))
    cat("Producing ", k, " random sets of variants and ", length(pval_cols)-k,
        " sets of variants using p-values in the data.\n")
    nm <- paste0(sample(c(letters, LETTERS), size=4, replace=TRUE), collapse="")
    nms <- paste0(nm, seq(k))
    pval_cols[is.na(pval_cols)] <- nms
    for(nm in nms) variants[[nm]] <- sample(seq(n), size=n, replace=FALSE)/n
  }
  stopifnot(all(pval_cols %in% names(variants)))

  variants <- variants %>% rename(snp = variant_name) %>%
              filter(snp %in% total_ld_variants )
  n <- nrow(variants)
  cat("Of these, ", n , " have LD information.\n")
  if(nrow(variants) == 0) return(NULL)

  ids <- variants$snp

  ld <- ld %>% filter(colsnp %in% ids & rowsnp %in% ids) %>%
        filter(r2 >= r2_thresh)
  o_c <- lapply(seq_along(pval_cols), function(i){
    nm <- pval_cols[i]
    ct <- sum(variants[[nm]]  < pval_thresh[i])
    order(variants[[nm]], decreasing=FALSE)[seq(ct)]
  })
  k_c <- lapply(seq_along(pval_cols), function(x){c()})
  for(j in seq_along(pval_cols)){
    while(length(o_c[[j]]) > 0){
      snp <- ids[o_c[[j]][1]]
      k_c[[j]] <- c(k_c[[j]], snp)
      if(length(o_c[[j]]) == 1){
        o_c[[j]] <- numeric()
        next
      }
      myld <- filter(ld, rowsnp==snp | colsnp==snp)
      if(nrow(myld) > 0){
        remove_snps <- with(myld, unique( c(rowsnp,colsnp)))
        remove_ix <- which(ids %in% remove_snps)
      }else{
        remove_ix <- o_c[[j]][1]
      }
      o_c[[j]] <- o_c[[j]][!o_c[[j]] %in% remove_ix]
      #cat(length(o_c), " ", length(o_c[[j]]), "\n")
    }
  }
  if(length(pval_cols)==1) return(unlist(k_c))
  return(k_c)
}
