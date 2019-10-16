library(dplyr)
library(readr)
library(MRPRESSO)

args <- commandArgs(trailingOnly=TRUE)
dat <- readRDS(args[1])
snp_list <- read_lines(args[2])
pval <- as.numeric(args[3])
out <- args[4]

dat <- filter(dat, p_value < pval & 
                   snp %in% snp_list)
dat <- data.frame(dat)
nbd <- max(1000, nrow(dat)*10)

#if(nrow(dat) < 3){
#    f1 <- NULL
#}else{
    f1 <- try(mr_presso(BetaOutcome = "beta_hat_2", BetaExposure = "beta_hat_1", 
              SdOutcome = "seb2", SdExposure = "seb1", OUTLIERtest = TRUE, 
              DISTORTIONtest = TRUE, data = dat, 
              NbDistribution = nbd,  SignifThreshold = 0.05), silent=TRUE)
#}
saveRDS(f1, file=out)
