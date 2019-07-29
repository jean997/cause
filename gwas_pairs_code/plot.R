library(tidyverse)
library(gridExtra)
library(cause)


args <- commandArgs(trailingOnly=TRUE)
cause_dir <- args[1]
mr_dir <- args[2]


tags <- c("glg_tg", "glg_ldl", "glg_hdl", "glg_tc",
          "giant_height", "giant_bmi", "lu_bfp",
          "egg_bl", "egg_bw", "egg_hc",
          #"icbp_dbp", "icbp_sbp", "icbp_map", "icbp_pp",
          "vanderHarst_cad", "megastroke_as",
          "gefos_bone", "ckdgen_egfrcrea", "magic_fg", "diagram_t2d")

trait_names = c("trig", "ldl", "hdl", "total chol",
                "height", "bmi", "body fat %",
                "birth ln", "birth wt",
                "head circ",
                #"dbp", "sbp", "map", "pulse press",
                "cad",  "stroke",
                "bone dens",
                "egfr crea",
                "fasting glucose", "t2d" )



tag_names <- data.frame(tag = tags,
                        name = trait_names)

tag_pairs <- expand.grid(tags, tags) %>%
            filter(Var1 != Var2) %>%
            mutate(fn_cause = paste0(cause_dir, Var1, "__", Var2, "_cause.RDS"),
                   fn_mr = paste0(mr_dir, Var1, "__", Var2, "_ivw.RDS"),
                   fn_mregger = paste0(mr_dir, Var1, "__", Var2, "_mregger.RDS"),
                   fn_mrp = paste0(mr_dir, Var1, "__", Var2, "_mrpresso.RDS")) %>%
            rename(tag1 = Var1,
                   tag2 = Var2)

tag_pairs$trait1 <- tag_names$name[match(tag_pairs$tag1, tag_names$tag)]
tag_pairs$trait2 <- tag_names$name[match(tag_pairs$tag2, tag_names$tag)]

tag_pairs$trait1 <- factor(tag_pairs$trait1, levels=trait_names)
tag_pairs$trait2 <- factor(tag_pairs$trait2, levels=trait_names)

dfcause <-data.frame(t(apply(tag_pairs, 1,  function(x){
             if(!file.exists(x["fn_cause"])){
                xx <- rep(NA, 4)
                names(xx) <- c("z_cause", "gamma_med3", "q_med2", "eta_med2")
                return(xx)
             }
             #cat(x["fn_cause"], "\n")
             res <- readRDS(x["fn_cause"])
             res_sum <- summary(res)
             z <- res_sum$z

             gamma_median3 <- res_sum$quants[[2]][1,1]
             q_median2 <- res_sum$quants[[1]][1,3]
             eta_median2 <- res_sum$quants[[1]][1,2]

             xx <- c(z, gamma_median3, q_median2, eta_median2)
             names(xx) <- c("z_cause", "gamma_med3", "q_med2", "eta_med2")
            return(xx)
       })))


df_mr <-data.frame(t(apply(tag_pairs, 1,  function(x){
             ivw <- readRDS(x["fn_mr"])
             egger <- readRDS(x["fn_mregger"])
             mrp <- readRDS(x["fn_mrp"])
             if(is.null(ivw)){
                 p_ivw <- NA
                 est_ivw <- NA
             }else{
                 p_ivw <- as.numeric(summary(ivw)$coefficients[1,4])
                 est_ivw <- as.numeric(summary(ivw)$coefficients[1,1])
             }
             if(is.null(egger)){
                 p_egger <- NA
                 est_egger <- NA
             }else{
                 p_egger <- as.numeric(summary(egger)$coefficients[2,4])
                 est_egger <- as.numeric(summary(egger)$coefficients[2,1])
             }
             if(is.null(mrp) | class(mrp) == "try-error"){
                 p_mrp <- NA
                 est_mrp <- NA
             }else{
                r <- mrp["Main MR results"]$`Main MR results`[2, c(3,6)]
                if(is.na(r[1])){
                    r <- mrp["Main MR results"]$`Main MR results`[1, c(3,6)]
                }
                p_mrp <- as.numeric(r[2])
                est_mrp <- as.numeric(r[1])
             }
             xx <- c(p_ivw,est_ivw,  p_egger, est_egger,  p_mrp, est_mrp)
             names(xx) <- c("p_ivw", "est_ivw",  "p_egger", "est_egger", "p_mrp", "est_mrp")
            return(xx)
            })))

dfcause <- dfcause %>% mutate(p_cause = pnorm(z_cause, lower.tail=TRUE))
df <- cbind(tag_pairs, dfcause, df_mr)



df <- mutate(df, lab_q = case_when(q_med2 < qbeta(0.5, 1, 10) ~ "",
                                     TRUE ~ as.character(round(q_med2, digits=2))),
                 fdr_cause = p.adjust(df$p_cause, method="BH"),
                 fdr_ivw = p.adjust(df$p_ivw, method="BH"),
                 fdr_egger = p.adjust(df$p_egger,  method="BH"),
                 fdr_mrp = p.adjust(df$p_mrp,  method="BH"))

saveRDS(df, file=paste0(cause_dir, "df.RDS"))
col1 <- "darkorange2"
col2 <- "cornflowerblue"
#### Combo plot ####


X1 <- select(df, fdr_cause, trait1, trait2, gamma_med3) %>%
        mutate(sign = sign(gamma_med3) , t1 = paste0(trait1, "_cause"), method="cause") %>%
        select(-gamma_med3, -trait1) %>%
        rename(fdr = fdr_cause, trait1 = t1) %>%
        select(trait1, trait2, fdr, sign, method)

X2 <- select(df, fdr_ivw, trait1, trait2, est_ivw) %>%
        mutate(sign = sign(est_ivw) , t1 = paste0(trait1, "_mr"), method="ivw") %>%
        select(-est_ivw, -trait1) %>%
        rename(fdr = fdr_ivw, trait1 = t1) %>%
        select(trait1, trait2, fdr, sign, method)

X <- rbind(X1, X2)

t1_order <- paste0(rep(trait_names, each=2), c("_cause", "_mr"))
X$trait1 <- factor(X$trait1, levels=t1_order)
plt_combo <- filter(X, fdr < 0.05) %>%
             ggplot(.) +
                geom_point(aes(x=as.numeric(trait2), y=as.numeric(trait1),
                                            size = -log10(fdr), fill=factor(method),
                                            color=factor(method), shape = factor(sign)))  +
                xlab("Outcome") + ylab("Mediator") +
                scale_color_manual(values=c("cause" = col2, "ivw" =col1), name="Method", labels=c("CAUSE", "IVW")) +
                scale_fill_manual(values=c("cause" = col2, "ivw" =col1), name="Method", labels=c("CAUSE", "IVW")) +
                scale_shape_manual(values=c("-1" = 25, "1" = 24), name="Effect Direction", labels=c("-", "+")) +
                scale_size_continuous(name="-log10(q-value)") +
                scale_y_continuous(breaks = seq(0.5, (length(trait_names)*2 - 1.5), by=2),
                                   labels=trait_names, limits=c(0.5, 2*length(trait_names)+0.5),
                                   expand = c(0, 0)) +
                scale_x_continuous(breaks = seq(0.5, length(trait_names)-0.5, by=1),
                                   labels=trait_names, limits=c(0.5, length(trait_names)+0.5),
                                   expand=c(0,0)) +
                geom_point(aes(x=as.numeric(trait2), y=as.numeric(trait1)), col="grey",
                           fill = "grey", shape=15, data = filter(X, is.na(fdr)), size=2) +
                geom_point(aes(x=x, y=y), shape=4, size=4,
                           data = data.frame(y = seq(1.5, 39.5, by=2), x = 1:20)) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 2),
                       panel.grid.minor=element_blank(),
                       axis.text.y = element_text(vjust = -1.5),
                       axis.ticks.y = element_blank())

ggsave(plt_combo, file=paste0(cause_dir, "gwas_fdr_cause_ivw.png"), height=8.5, width=8.5, units="in", dpi=300)

####################
# q plot

plt <-  filter(df, q_med2 > qbeta(0.5, 1, 10)) %>%
        ggplot(.) +
            geom_point(aes(x=as.numeric(trait2), y=as.numeric(trait1), shape=factor(sign(eta_med2), levels=c("-1", "1")),
                               fill=factor(fdr_cause < 0.05, levels=c(TRUE, FALSE)),
                               color=factor(fdr_cause < 0.05, levels=c(TRUE, FALSE)),
                               size=q_med2)) +
            geom_text( aes(x=as.numeric(trait2), y=as.numeric(trait1) + 0.3, label=lab_q), size=2) +
            xlab("Outcome") + ylab("Mediator") +
            scale_shape_manual(values=c("-1"=25, "1"=24), name="Effect Direction", labels=c("-", "+")) +
            scale_fill_manual(values=c("FALSE"="white", "TRUE"="grey"), name="Sharing Model\nRejected", guide=FALSE) +
            scale_color_manual(values=c("FALSE"="black", "TRUE"="white"), name="Sharing Model\nRejected", guide=FALSE) +
            scale_size_continuous(name="Proportion Shared Variants") +
            scale_y_continuous(breaks = seq(0.5, length(trait_names)-0.5, by=1),
                                   labels=trait_names, limits=c(0.5, length(trait_names)+0.5),
                                   expand=c(0,0)) +
            scale_x_continuous(breaks = seq(0.5, length(trait_names)-0.5, by=1),
                                   labels=trait_names, limits=c(0.5, length(trait_names)+0.5),
                                   expand=c(0,0)) +
            geom_point(aes(x=x, y=y), shape=4, size=4,
                           data = data.frame(y = 1:20, x = 1:20)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 2),
                  panel.grid.minor=element_blank(),
                  axis.text.y = element_text(vjust = -1.5),
                  axis.ticks.y = element_blank())

ggsave(plt, file=paste0(cause_dir, "gwas_q_prop_m2.png"), height=8, width=10, units="in", dpi=300)


############## All MR Methods together


X1 <- select(df, fdr_mrp, trait1, trait2, est_mrp) %>%
        mutate(sign = sign(est_mrp) , t1 = paste0(trait1, "_mrp"), method="MR-PRESSO") %>%
        select(-est_mrp, -trait1) %>%
        rename(fdr = fdr_mrp, trait1 = t1) %>%
        select(trait1, trait2, fdr, sign, method)

X2 <- select(df, fdr_ivw, trait1, trait2, est_ivw) %>%
        mutate(sign = sign(est_ivw) , t1 = paste0(trait1, "_mr"), method="IVW") %>%
        select(-est_ivw, -trait1) %>%
        rename(fdr = fdr_ivw, trait1 = t1) %>%
        select(trait1, trait2, fdr, sign, method)

X3 <- select(df, fdr_egger, trait1, trait2, est_egger) %>%
        mutate(sign = sign(est_egger) , t1 = paste0(trait1, "_egger"), method="Egger") %>%
        select(-est_egger, -trait1) %>%
        rename(fdr = fdr_egger, trait1 = t1) %>%
        select(trait1, trait2, fdr, sign, method)


X <- rbind(X1, X2, X3)

t1_order <- paste0(rep(trait_names, each=3), c("_mrp", "_egger", "_mr"))
X$trait1 <- factor(X$trait1, levels=t1_order)
X$method <- factor(X$method, levels=c("IVW", "Egger", "MR-PRESSO"))

#col3 <- "chartreuse3"
col3 <- "gold"
col4 <- "darkorchid"

plt_combo <- filter(X, fdr < 0.05 & !is.na(fdr)) %>%
             ggplot(.) +
                geom_point(aes(x=as.numeric(trait2), y=as.numeric(trait1),
                               size = -log10(fdr), fill=factor(method),
                               color=factor(method), shape = factor(sign))) +
                xlab("Outcome") + ylab("Mediator") +
                scale_color_manual(values=c("MR-PRESSO" = col4, "Egger"=col3, "IVW" =col1), name="Method") +
                scale_fill_manual(values=c("MR-PRESSO" = col4, "Egger"=col3, "IVW" =col1), name="Method") +
                scale_shape_manual(values=c("-1" = 25, "1" = 24), name="Effect Direction", labels=c("-", "+")) +
                scale_size_continuous(name="-log10(q-value)") +
                scale_y_continuous(breaks = seq(0.5, (length(trait_names)*3 - 1.5), by=3),
                                   labels=trait_names, limits=c(0.5, 3*length(trait_names) + 0.5),
                                   expand=c(0, 0)) +
                scale_x_continuous(breaks = seq(0.5, length(trait_names)-0.5, by=1), labels=trait_names,
                                   limits=c(0.5, length(trait_names) + 0.5), expand=c(0, 0)) +
                geom_point(aes(x=as.numeric(trait2), y=as.numeric(trait1)), col="grey",
                           fill = "grey", shape=15, data = filter(X, is.na(fdr)), size=2) +
                geom_point(aes(x=x, y=y), shape=4, size=4,
                           data = data.frame(y = seq(2, 59, by=3), x = 1:20)) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 2),
                       panel.grid.minor=element_blank(),
                       axis.text.y = element_text(vjust = -1.5),
                       axis.ticks.y = element_blank())

ggsave(plt_combo, file=paste0(cause_dir, "gwas_fdr_ivw_egger_mrp.png"), height=10, width=8.5, units="in", dpi=300)

