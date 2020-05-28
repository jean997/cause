library(tidyverse)
library(cause)


args <- commandArgs(trailingOnly=TRUE)
#Input files
data_file <- args[1]
snp_file_asc <- args[2]
params <- args[3]
# q prior parameters
qalpha <- as.numeric(args[4])
qbeta <- as.numeric(args[5])
max_q <- as.numeric(args[6])
#Output files
cause_out <- args[7]
seed <- as.numeric(args[8])

X <- readRDS(data_file)
params <- readRDS(params)
snps_asc <- read_lines(snp_file_asc)

set.seed(seed)
res <- cause(X =X, param_ests = params, variants=snps_asc,
             qalpha = qalpha, qbeta = qbeta, max_q = max_q,  force=TRUE)
saveRDS(res, cause_out)

