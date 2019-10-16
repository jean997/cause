library(dplyr)
library(readr)
library(MendelianRandomization)


args <- commandArgs(trailingOnly=TRUE)
file1 <- args[1]
file2 <- args[2]
snp_file_asc <- args[3]
pval_thresh <- as.numeric(args[4])
output_file <- args[5]

X1 <- cause::read_standard_format(file1)
X2 <- cause::read_standard_format(file2)

X <- cause::gwas_merge(X1, X2, X1_formatted=TRUE, X2_formatted = TRUE)
X$p_value <- with(X1, p_value[match(X$snp, snp)])

snp_list <- read_lines(snp_file_asc)

results <- tibble(method=c("IVW_RE_noNOME", "Egger_RE_NOME", 
                           "Median_Wtd", "MBE_Wtd_noNOME_phi1",
                           "MBE_Wtd_noNOME_phi0.5", 
                           "MBE_Wtd_noNOME_phi0.25"), 
                  est = NA, se = NA, pvalue = NA) 
if(with(X, sum(p_value < pval_thresh & snp %in% snp_list) >= 3)){
    dat <- filter(X, p_value < pval_thresh & 
                   snp %in% snp_list) %>%
          with(., mr_input(bx = beta_hat_1, bxse = seb1, 
                            by = beta_hat_2, byse = seb2, 
                            snps = snp))


    f <- mr_ivw(dat, model="random", weights="delta")
    results$est[1] <- f@Estimate
    results$se[1] <- f@StdError
    results$pvalue[1] <- f@Pvalue


    f <- mr_egger(dat)
    results$est[2] <- f@Estimate
    results$se[2] <- f@StdError.Est
    results$pvalue[2] <- f@Pvalue.Est

    f <- mr_median(dat)
    results$est[3] <- f@Estimate
    results$se[3] <- f@StdError
    results$pvalue[3] <- f@Pvalue

    f <- mr_mbe(dat, stderror = "delta", phi=1)
    results$est[4] <- f@Estimate
    results$se[4] <- f@StdError
    results$pvalue[4] <- f@Pvalue


    f <- mr_mbe(dat, stderror = "delta", phi=0.5)
    results$est[5] <- f@Estimate
    results$se[5] <- f@StdError
    results$pvalue[5] <- f@Pvalue


    f <- mr_mbe(dat, stderror = "delta", phi=0.25)
    results$est[6] <- f@Estimate
    results$se[6] <- f@StdError
    results$pvalue[6] <- f@Pvalue
}
saveRDS(results, file=output_file)

