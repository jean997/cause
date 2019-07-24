library(tidyverse)

res_file <- "sigma_g_causedf.RDS"
causedf <- readRDS(res_file)

vars <- c("z", "p", "eta_sharing", "q_sharing", 
          "eta_causal", "q_causal", "gamma_causal")

min_cor <- expand.grid(variable = vars, 
                       params = unique(causedf$params), stringsAsFactors=F)
min_cor$min_cor <- NA
min_cor$cor_med_small <- NA
for(i in seq(nrow(min_cor))){
    v <- min_cor$variable[i]
    X <- select(causedf, simulate.output.file, params, quant, one_of(v)) %>%
         spread(quant, get(v))  %>%
         filter(params==min_cor$params[i]) %>%
         select(one_of(c("0.51", "0.65", "0.8")))

    min_cor$min_cor[i] <- min(cor(X))
    min_cor$cor_med_small[i] <- cor(X[,"0.65"], X[,"0.8"])

}


