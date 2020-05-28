library(dplyr)
library(readr)
library(cause)


args <- commandArgs(trailingOnly=TRUE)
data <- args[1]
ld_score_dir <- args[2]
output_file <- args[3]

X <- readRDS(data)

X_ldsc <- purrr::map_df(1:22, function(i){
    cat(i, " ")
    ldscores <- read_tsv(paste0(ld_score_dir, i, ".l2.ldscore.gz")) %>%
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
