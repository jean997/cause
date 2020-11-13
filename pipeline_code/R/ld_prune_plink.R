library(tidyverse)
library(ieugwasr)

args <- commandArgs(trailingOnly=TRUE)
X <- readRDS(args[1])
pval_thresh <- as.numeric(args[2])
r2_thresh <- as.numeric(args[3])
ref_path <- args[4]
out_file <- args[5]

X <- X %>% rename(rsid = snp,
                  pval = p_value)

X_clump <- ld_clump(dat = X,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = ref_path)
keep_snps <- X_clump$rsid
saveRDS(keep_snps, file=out_file)
