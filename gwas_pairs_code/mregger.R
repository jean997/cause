library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
dat <- readRDS(args[1])
pval <- as.numeric(args[2])
out <- args[3]

dat <- filter(dat, -abs(beta_hat_1/seb1) < qnorm(pval/2) & 
                   seb1 > 0 & seb2 > 0)
if(nrow(dat) < 3){
    f1 <- NULL
}else{
    f1 <- lm(beta_hat_2 ~beta_hat_1, weights=1/(seb2^2), data=dat)
}
saveRDS(f1, file=out)
