library(dplyr)
library(readr)
library(cause)


args <- commandArgs(trailingOnly=TRUE)
file1 <- args[1]
file2 <- args[2]
output_file <- args[3]

X1 <- read_standard_format(file1)
X2 <- read_standard_format(file2)

X <- gwas_merge(X1, X2, X1_formatted=TRUE, X2_formatted = TRUE)
X$p_value <- with(X1, p_value[match(X$snp, snp)])

#hm3_snps <- read_tsv("/project2/xinhe/ldsc_reference_files/eur_w_ld_chr/w_hm3.snplist")

X_ldsc <- purrr::map_df(1:22, function(i){
    cat(i, " ")
    ldscores <- read_tsv(paste0("/project2/xinhe/ldsc_reference_files/eur_w_ld_chr/", i, ".l2.ldscore.gz")) %>%
                select(SNP, L2, CHR, BP) %>%
                rename(snp = SNP)
    inner_join(X, ldscores, by="snp")
})
o <- with(X_ldsc, order(CHR, BP))
X_ldsc <- X_ldsc[o,]
lcv_res <- with(X_ldsc, 
                causeSims::RunLCV(L2,
                       beta_hat_1/seb1,beta_hat_2/seb2,
                       no.blocks=100,crosstrait.intercept=1,ldsc.intercept=1,
                       sig.threshold=30,n.1=1,n.2=1,intercept.12=0))

saveRDS(lcv_res, file=output_file)
