library(tidyverse)
library(cause)


args <- commandArgs(trailingOnly=TRUE)
#Input files
file1 <- args[1]
file2 <- args[2]
snp_file_asc <- args[3]
params <- args[4]
# q prior parameters
qalpha <- as.numeric(args[5])
qbeta <- as.numeric(args[6])
max_q <- as.numeric(args[7])
#Output files
cause_out <- args[8]
seed <- as.numeric(args[9])

X1 <- read_standard_format(file1)
X2 <- read_standard_format(file2)

X <- gwas_merge(X1, X2, X1_formatted=TRUE, X2_formatted = TRUE)
X$p_value <- with(X1, p_value[match(X$snp, snp)])
#X <- readRDS(data_file)
params <- readRDS(params)

snps_asc <- read_lines(snp_file_asc)
res <- cause(X =X, param_ests = params, variants=snps_asc,
             qalpha = qalpha, qbeta = qbeta, max_q = max_q,  force=TRUE)
saveRDS(res, cause_out)

