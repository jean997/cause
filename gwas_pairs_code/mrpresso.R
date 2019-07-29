library(dplyr)
library(MRPRESSO)

args <- commandArgs(trailingOnly=TRUE)
dat <- readRDS(args[1])
pval <- as.numeric(args[2])
out <- args[3]

dat <- filter(dat, -abs(beta_hat_1/seb1) < qnorm(pval/2) & 
                   seb1 > 0 & seb2 > 0)
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
