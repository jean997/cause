library(tidyverse)
library(cause)

args <- commandArgs(trailingOnly=TRUE)
X <- readRDS(args[1])
pval_thresh <- as.numeric(args[2])
r2_thresh <- as.numeric(args[3])
ld <- readRDS(args[4])
snp_info <- readRDS(args[5])
out_file <- args[6]

keep_snps <- ld_prune(variants = X, ld, total_ld_variants = snp_info$SNP, variant_name="snp",
                  r2_thresh = r2_thresh, pval_cols = c("p_value"), pval_thresh = c(pval_thresh))

saveRDS(keep_snps, file=out_file)
