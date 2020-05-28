library(tidyverse)
library(cause)


args <- commandArgs(trailingOnly=TRUE)
#Input files
data_file <- args[1]
#Output file
params_out <- args[2]
#Seed
seed <- as.numeric(args[3])


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

