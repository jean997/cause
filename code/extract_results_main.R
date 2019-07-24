library(tidyverse)
library(dscrutils)


# Directories housing our results
dirs <- c("pwr", "fp")
labs <- c("pwr", "fp")

# For each directory, extract a data frame of simulation parameters
# params is a list of two data frames
params <- lapply(1:2, FUN = function(i){
    dir <- dirs[i]
    lab <- labs[i]
    paramdf <- dscquery(dsc.outdir=dir, 
                   targets=c("simulate", 
                             "gw_sig.tau","gw_sig.gamma", "gw_sig.q", 
                             "gw_sig.omega", "gw_sig.eta", "gw_sig.h1", "gw_sig.h2",
                             "gw_sig.n1", "gw_sig.n2", "gw_sig.neffect1",
                             "gw_sig.neffect2", "gw_sig.lcv_q1", "gw_sig.lcv_q2",
                             "gw_sig.lcv_gcp",  
                             "gw_sig.m_sig") , omit.filenames=F)
     # aesthetic re-naming                         
     paramdf <- paramdf %>% rename(q = gw_sig.q, tau = gw_sig.tau, 
                                  eta = gw_sig.eta, gamma = gw_sig.gamma,
                                  omega = gw_sig.omega, h1 = gw_sig.h1, 
                              h2 = gw_sig.h2, n1 = gw_sig.n1, n2 = gw_sig.n2, 
                              neffect1 = gw_sig.neffect1, neffect2 = gw_sig.neffect2,
                              m_sig = gw_sig.m_sig, 
                              lcv_q1 = gw_sig.lcv_q1, lcv_q2 = gw_sig.lcv_q2, 
                              gcp = gw_sig.lcv_gcp) %>%
                mutate(params = paste0("(", tau, ",", omega, ",", q, ")")) 
     return(paramdf)
    })

# For each directory extract a data frame of CAUSE results
# cause_res is a list of two data frames
# We will extract p-values but also posterior medians from sharing and causal models
cause_res <- lapply(1:2, FUN = function(i){
    dir <- dirs[i]
    lab <- labs[i]
    causedf <- dscquery(dsc.outdir=dir, 
                   targets=c("simulate", "cause",
                             "cause.p", "cause.eta_med_2", "cause.q_med_2", 
                             "cause.eta_med_3", "cause.q_med_3", "cause.gamma_med_3"),
                      ignore.missing.files=T,
                      omit.filenames=F) 

    causedf <- causedf %>% 
               rename(cause.eta_sharing = cause.eta_med_2, cause.q_sharing = cause.q_med_2,
                       cause.eta_causal = cause.eta_med_3, cause.q_causal = cause.q_med_3,
                       cause.gamma_causal = cause.gamma_med_3)  %>%
             full_join(params[[i]], ., by=c("DSC", "simulate.output.file")) %>%
             mutate(tag = paste0(lab, "_", simulate.output.file))
    return(causedf)
    })
# Combine the two data frames into one and save it
causedf <- bind_rows(cause_res[[1]], cause_res[[2]])
saveRDS(causedf, paste0("main_causedf_", Sys.Date(), ".RDS")) 

# For each directory extract a data frame of LCV results
# cause_res is a list of two data frames
# We extract p-values and the posterior mean and standard error of the GCP
lcv_res <- lapply(1:2, FUN = function(i){
    dir <- dirs[i]
    lab <- labs[i]
    lcvdf <- dscquery(dsc.outdir=dir, 
         targets=c("simulate", 
                   "LCV.p", "LCV.gcp_mean", "LCV.gcp_pse"), 
                   omit.filenames=F)

    lcvdf <- lcvdf %>% 
             full_join(params[[i]], ., by=c("DSC", "simulate.output.file")) %>%
             mutate(tag = paste0(lab, "_", simulate.output.file))
    return(lcvdf)
    })
# Combine the two data frames into one and save it
lcvdf <- bind_rows(lcv_res[[1]], lcv_res[[2]])
saveRDS(lcvdf, paste0("main_lcvdf_", Sys.Date(), ".RDS")) 

# For each directory extract a data frame of other MR results
# mr_res is a list of two data frames
# For these methods we only extrac the p-value
mr_res <- lapply(1:2, FUN = function(i){
    dir <- dirs[i]
    lab <- labs[i]
    mrdf <- dscquery(dsc.outdir=dir, 
                   targets=c("simulate", 
                             "mr.p"),
                      omit.filenames=F)

    mrdf <- mrdf %>% 
             full_join(params[[i]], ., by=c("DSC", "simulate.output.file")) %>%
             mutate(tag = paste0(lab, "_", simulate.output.file))
    return(mrdf)
    })
mrdf <- bind_rows(mr_res[[1]], mr_res[[2]])
saveRDS(mrdf, paste0("main_mrdf_", Sys.Date(), ".RDS")) 



