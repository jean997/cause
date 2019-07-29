library(tidyverse)
library(cause)

args <- commandArgs(trailingOnly=TRUE)
data1 <- args[1]
chrom <- as.numeric(args[2])
pval_thresh <- as.numeric(args[3])
r2_thresh <- as.numeric(args[4])
ld <- readRDS(args[5])
snp_info <- readRDS(args[6])
out_file <- args[7]

d1 <- read_tsv(data1, col_types = list(col_character(), col_integer(),
                                       col_character(),col_character(),
                                       col_character(), col_double(),
                                       col_double(),col_double(),
                                       col_character(),col_double())) 

dup1 <- d1$snp[duplicated(d1$snp)]
d1 <- d1 %>% filter(!snp %in% dup1)

keep_snps <- ld_prune(variants = d1, ld, total_ld_variants = snp_info$SNP, variant_name="snp",
                  r2_thresh = r2_thresh, pval_cols = c("p_value"), pval_thresh = c(pval_thresh))

saveRDS(keep_snps, file=out_file)
