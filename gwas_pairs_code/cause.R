library(tidyverse)
library(cause)


args <- commandArgs(trailingOnly=TRUE)
#Input files
data_file_1 <- args[1]
data_file_2 <- args[2]
snp_file_asc <- args[3]
#Output files
params_out <- args[4]
cause_out <- args[5]
data_out <- args[6]
seed <- as.numeric(args[7])

#if(is.na(seed)) seed <- 100


d1 <- read_tsv(data_file_1, col_type=list(col_character(),col_integer(),
                                 col_character(),col_character(),
                                 col_character(),col_double(),col_double(),
                                 col_double(),col_double(),
                                 col_double())) 

dup1 <- d1$snp[duplicated(d1$snp)]
d1 <- d1 %>% filter(!snp %in% dup1)

d2 <- read_tsv(data_file_2, col_type=list(col_character(),col_integer(),
                                 col_character(),col_character(),
                                 col_character(),col_double(),col_double(),
                                 col_double(),col_double(),
                                 col_double())) 

dup2 <- d2$snp[duplicated(d2$snp)]
d2 <- d2 %>% filter(!snp %in% dup2)

#This what we would use if the data weren't already aligned and strand flipped uniformly
# It takes about 9 seconds
#X <- gwas_format_cause(d1, d2, snp_name_cols = rep("snp", 2),
#                       beta_hat_cols = rep("beta_hat", 2),
#                       se_cols = rep("se", 2), A1_cols = rep("ref_allele", 2),
#                       A2_cols = rep("alt_allele", 2))

# Instead we can use this since the data are cleaned up uniformly
# This takes about 1 second because it doesn't have to check for strand flips
X <- d1 %>%
     select(snp, beta_hat, se) %>%
     rename(beta_hat_1 = beta_hat, seb1 = se) %>%
     inner_join(., d2, by="snp") %>%
     rename(beta_hat_2 = beta_hat, seb2 = se,
           A1 = ref_allele, A2 = alt_allele) %>%
     select(snp, beta_hat_1, seb1, beta_hat_2, seb2, A1, A2) %>%
     new_cause_data(.)

# We save the processed data for use with other MR methods
saveRDS(X, file=data_out)

set.seed(seed)
if(nrow(X) < 1e6){
    snps_grid <- X$snp
}else{
    snps_grid <- sample(X$snp, size=1e6, replace=FALSE)
}

t <- system.time(params <- est_cause_params(X, snps_grid))
params$time <- t
saveRDS(params, params_out)

snps_asc <- read_lines(snp_file_asc)
res <- cause(X =X, param_ests = params, variants=snps_asc,
             qalpha = 1, qbeta = 10, force=TRUE)
saveRDS(res, cause_out)

