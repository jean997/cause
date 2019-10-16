library(stringr)
library(dplyr)
library(purrr)
library(cause)

args <- commandArgs(trailingOnly=TRUE)
results_dir <- args[1]
method <- args[2]

fl <- list.files(results_dir, paste0(method, ".RDS$"), full.names=TRUE)
trait_pairs <- str_replace(fl, paste0(results_dir, "/"), "")
trait_pairs <- str_replace(trait_pairs,  paste0("_", method, ".RDS$"), "")
trait_pairs <- str_split(trait_pairs, "__")
tag1 <- purrr::map(trait_pairs, 1) %>% unlist()
tag2 <- purrr::map(trait_pairs, 2) %>% unlist()


df  <- tibble(tag1 = tag1, tag2 = tag2, file=fl)

if(method=="ivw"){
    dfres <-map_dfr(df$file, function(f){
             res <- readRDS(f)
             if(is.null(res)){
                 p <- NA
                 est <- NA
             }else{
                 p <- as.numeric(summary(res)$coefficients[1,4])
                 est <- as.numeric(summary(res)$coefficients[1,1])
             }
             xx <- tibble(file=f,  estimate=est, p=p)
            return(xx)
       }) 
       dfres$method <- method
}else if(method=="egger"){
    dfres <-map_dfr(df$file, function(f){
             res <- readRDS(f)
             if(is.null(res)){
                 p <- NA
                 est <- NA
             }else{
                 p <- as.numeric(summary(res)$coefficients[2,4])
                 est <- as.numeric(summary(res)$coefficients[2,1])
             }
             xx <- tibble(file=f,  estimate=est, p=p)
            return(xx)
       })
       dfres$method <- method
}else if(method=="mrpresso"){
    dfres <-map_dfr(df$file, function(f){
             #cat(f, "\n")
             res <- readRDS(f)
             if(is.null(res) | class(res) == "try-error"){
                 p <- NA
                 est <- NA
             }else{
                r <- res["Main MR results"]$`Main MR results`[2, c(3,6)]
                if(is.na(r[1])){
                    r <- res["Main MR results"]$`Main MR results`[1, c(3,6)]
                }
                p <- as.numeric(r[2])
                est <- as.numeric(r[1])
             }
             xx <- tibble(file=f,  estimate=est, p=p)
            return(xx)
       })
       dfres$method <- method
}else if(method=="mbe"){
    dfres <- map_dfr(df$file, function(f){
                           res <- readRDS(f)
                           if(is.null(res)){
                               return(NULL)
                            }
                            res <- res %>% mutate(file = f, 
                                                  method = paste0("mbe_", Method, "_", phi)) %>%
                                    select(-Method, -phi, -CI_low, -CI_upp)

                            return(res)
                        }) 
}else if(str_starts(method, "cause")){
    dfres <-map_dfr(df$file, function(f){
             res <- readRDS(f)
             res_sum <- summary(res)
             z <- res_sum$z

             gamma_median3 <- res_sum$quants[[2]][1,1]
             q_median3 <- res_sum$quants[[2]][1,3]
             eta_median3 <- res_sum$quants[[2]][1,2]
             q_median2 <- res_sum$quants[[1]][1,3]
             eta_median2 <- res_sum$quants[[1]][1,2]

             xx <- tibble(file=f, z=z, gamma_med3=gamma_median3, 
                          q_me3 = q_median3, eta_med3 = eta_median3,
                          q_med2 = q_median2, eta_med2 = eta_median2)
            return(xx)
       })
    dfres$method <- method
}else if(method=="lcv"){
    dfres <-map_dfr(df$file, function(f){
             res <- readRDS(f)
             if(is.null(res)){
                 p <- NA
                 z <- NA
                 gcp_pm <- NA
                 gcp_pse <- NA
                 rho <- rho_se <- NA
             }else{
                 p <- res$pval.gcpzero.2tailed
                 z <- res$zscore
                 gcp_pm <- res$gcp.pm
                 gcp_pse <- res$gcp.pse
                 rho <- res$rho.est
                 rho_se <- res$rho.err
             }
             xx <- tibble(file=f,  p=p, z=z, 
                          gcp_pm = gcp_pm, gcp_pse= gcp_pse, 
                          rho = rho, rho_se = rho_se)
            return(xx)
       })
    dfres$method <- method
}else if(method=="mrpackage"){
    dfres <- map_dfr(df$file, function(f){
                           res <- readRDS(f)
                           if(is.null(res)){
                               return(NULL)
                            }
                            res <- res %>% mutate(file = f )
                            return(res)
                        }) 


}
df <- left_join(df, dfres, by="file")

cat("saving", paste0(results_dir, "df_", method, ".RDS"), "\n")
saveRDS(df, file=paste0(results_dir, "df_", method, ".RDS"))
