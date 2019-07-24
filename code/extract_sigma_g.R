library(tidyverse)
library(dscrutils)


dir <- "sigma_g"
lab <- "sg"

paramdf <- dscquery(dsc.outdir=dir, 
                   targets=c("simulate", 
                             "gw_sig.tau","gw_sig.gamma", "gw_sig.q", 
                             "gw_sig.omega", "gw_sig.eta", "gw_sig.h1", "gw_sig.h2",
                             "gw_sig.n1", "gw_sig.n2", "gw_sig.neffect1",
                             "gw_sig.neffect2", "gw_sig.lcv_q1", "gw_sig.lcv_q2",
                             "gw_sig.lcv_gcp",  
                             "gw_sig.m_sig") , omit.filenames=F)
                            
paramdf <- paramdf %>% rename(q = gw_sig.q, tau = gw_sig.tau, 
                                  eta = gw_sig.eta, gamma = gw_sig.gamma,
                                  omega = gw_sig.omega, h1 = gw_sig.h1, 
                              h2 = gw_sig.h2, n1 = gw_sig.n1, n2 = gw_sig.n2, 
                              neffect1 = gw_sig.neffect1, neffect2 = gw_sig.neffect2,
                              m_sig = gw_sig.m_sig, 
                              lcv_q1 = gw_sig.lcv_q1, lcv_q2 = gw_sig.lcv_q2, 
                              gcp = gw_sig.lcv_gcp) %>%
                mutate(params = paste0("(", tau, ",", omega, ",", q, ")")) 


causedf <- dscquery(dsc.outdir=dir, 
                   targets=c("simulate", "cause_sigma_g",
                             "cause_sigma_g.quant", "cause_sigma_g.sigma_g",
                             "cause_sigma_g.z", "cause_sigma_g.eta_med_2", 
                             "cause_sigma_g.q_med_2", 
                             "cause_sigma_g.eta_med_3", "cause_sigma_g.q_med_3", 
                             "cause_sigma_g.gamma_med_3"),
                      omit.filenames=F)
causedf <- causedf %>% 
           full_join(paramdf, ., by=c("DSC", "simulate.output.file")) %>%
           mutate(cause_sigma_g.p = pnorm(cause_sigma_g.z, lower.tail=FALSE))


causedf <- causedf %>% 
           rename(quant = cause_sigma_g.quant, 
                  sigma_g = cause_sigma_g.sigma_g,
                  z = cause_sigma_g.z,
                  p = cause_sigma_g.p,
                  eta_sharing = cause_sigma_g.eta_med_2,
                  q_sharing = cause_sigma_g.q_med_2,
                  eta_causal = cause_sigma_g.eta_med_3,
                  q_causal = cause_sigma_g.q_med_3,
                  gamma_causal = cause_sigma_g.gamma_med_3)

saveRDS(causedf, paste0(lab, "_causedf_", Sys.Date(), ".RDS")) 

