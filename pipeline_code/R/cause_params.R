library(tidyverse)
library(cause)


args <- commandArgs(trailingOnly=TRUE)
#Input files
data_file <- args[1]
snp_file_asc <- args[2]
#Output file
params_out <- args[3]
#Seed
seed <- as.numeric(args[4])


X <- readRDS(data_file)

set.seed(seed)
if(nrow(X) < 1e6){
    snps_grid <- X$snp
}else{
    snps_grid <- sample(X$snp, size=1e6, replace=FALSE)
}

t <- system.time(params <- est_cause_params(X, snps_grid))
params$time <- t
saveRDS(params, params_out)

