library(dplyr)
library(cause)

args <- commandArgs(trailingOnly=TRUE)
file1 <- args[1]
file2 <- args[2]
output_file <- args[3]

X1 <- read_standard_format(file1)
X2 <- read_standard_format(file2)

X <- gwas_merge(X1, X2, X1_formatted=TRUE, X2_formatted = TRUE)
X$p_value <- with(X1, p_value[match(X$snp, snp)])
saveRDS(X, file=output_file)
