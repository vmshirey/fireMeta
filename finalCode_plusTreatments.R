# load required packages, establish re-utilized variables, custom functions
rm(list = ls())

require(tidyverse); require(bayesmeta); library(meta); library(ggmap); library(gridExtra); library(ggExtra); library(cowplot)
library(compute.es)
library(esc)

# drops a single effect from the analysis and re-runs the meta-analysis for sensitivity
jackknife <- function(effects){
  eff.copy <- data.frame()
  
  dropVec <- unique(effects[,"EffectCode"]) # publication labels
  sens <- list() # list for critical values
  for(i in 1:length(dropVec)){
    eff.copy <- effects %>% filter(EffectCode != dropVec[i])# drop random study
    
    # run meta on reduced dataset
    ma <- rma.mv(yi=EffectSize, 
                 V=EffectSizeStandardError,
                 random=~1|PublicationCode/EffectCode,
                 tdist=TRUE,
                 data=eff.copy, 
                 method="REML",
                 slab=EffectCode)
    
    # save variables of interest to list
    label <- dropVec[i]
    sens[[i]] <- c(ma$pval, ma$b, ma$se, ma$ci.lb, ma$ci.ub, ma$tau2, label)
  }
  
  return(sens)
}

plots.list <- list() # list for storing plots
scores.list <- list() # list for storing z-scores

dat <- read.csv("FinalEffects.csv", stringsAsFactors = FALSE, sep=",", header=TRUE)
dat <- dat %>% filter(FireType == "Treatment+Fire") 

ma <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=dat, 
               method="REML",
               slab=EffectCode)

funnel(ma, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

forest(ma, order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)

# subset data by taxon/diversity metric
babun <- dat %>% filter(Metric=="Abundance") # abundance
brich <- dat %>% filter(Metric=="Richness") # richness

ma <- rma.mv(yi=EffectSize, 
             V=EffectSizeStandardError,
             random=~1|PublicationCode/EffectCode,
             tdist=TRUE,
             data=brich, 
             method="REML",
             slab=EffectCode)

funnel(ma, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
forest(ma, order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)

ma <- rma.mv(yi=EffectSize, 
             V=EffectSizeStandardError,
             random=~1|PublicationCode/EffectCode,
             tdist=TRUE,
             data=babun, 
             method="REML",
             slab=EffectCode)

funnel(ma, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
forest(ma, order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)

# Abundance: Bees
ma01 <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=filter(babun, Taxon=="Bees"), 
               method="REML",
               slab=EffectCode)

# jackknife analysis for abundance
babun.sens1 <- jackknife(filter(babun, Taxon=="Bees"))
estimate1 <- as.numeric(sapply(babun.sens1, "[[", 2)) # effect estimate
estimate1.lb <- as.numeric(sapply(babun.sens1, "[[", 4)) # lower bound
estimate1.ub <- as.numeric(sapply(babun.sens1, "[[", 5)) # upper bound
pubs.out <- sapply(babun.sens1, "[[", 7) # publication

scores.list[[1]] <- (ma01$b - mean(estimate1))/sd(estimate1)

mean(estimate1.lb); mean(estimate1.ub)

# Richness: Bees
ma02 <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=filter(brich, Taxon=="Bees"), 
               method="REML",
               slab=EffectCode)

# jackknife analysis for abundance
brich.sens1 <- jackknife(filter(brich, Taxon=="Bees"))
estimate2 <- as.numeric(sapply(brich.sens1, "[[", 2)) # effect estimate
estimate2.lb <- as.numeric(sapply(brich.sens1, "[[", 4)) # lower bound
estimate2.ub <- as.numeric(sapply(brich.sens1, "[[", 5)) # upper bound
pubs.out <- sapply(brich.sens1, "[[", 7) # publication

scores.list[[2]] <- (ma02$b - mean(estimate2))/sd(estimate2)

mean(estimate2.lb); mean(estimate2.ub)

################## BUTTERFLIES ##############################
# Abundance: Butterflies
ma03 <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=filter(babun, Taxon=="Butterflies"), 
               method="REML",
               slab=EffectCode)

# jackknife analysis for abundance
babun.sens1 <- jackknife(filter(babun, Taxon=="Butterflies"))
estimate3 <- as.numeric(sapply(babun.sens1, "[[", 2)) # effect estimate
estimate3.lb <- as.numeric(sapply(babun.sens1, "[[", 4)) # lower bound
estimate3.ub <- as.numeric(sapply(babun.sens1, "[[", 5)) # upper bound
pubs.out <- sapply(babun.sens1, "[[", 7) # publication

scores.list[[3]] <- (ma03$b - mean(estimate3))/sd(estimate3)

mean(estimate3.lb); mean(estimate3.ub)

# Richness: Butterflies
ma04 <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=filter(brich, Taxon=="Butterflies"), 
               method="REML",
               slab=EffectCode)

# jackknife analysis for abundance
brich.sens1 <- jackknife(filter(brich, Taxon=="Butterflies"))
estimate4 <- as.numeric(sapply(brich.sens1, "[[", 2)) # effect estimate
estimate4.lb <- as.numeric(sapply(brich.sens1, "[[", 4)) # lower bound
estimate4.ub <- as.numeric(sapply(brich.sens1, "[[", 5)) # upper bound
pubs.out <- sapply(brich.sens1, "[[", 7) # publication

scores.list[[4]] <- (ma04$b - mean(estimate4))/sd(estimate4)

################## CARABIDAE ##############################
# Abundance: Carabids
ma05 <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=filter(babun, Taxon=="Carabids"), 
               method="REML",
               slab=EffectCode)

# jackknife analysis for abundance
babun.sens1 <- jackknife(filter(babun, Taxon=="Carabids"))
estimate5 <- as.numeric(sapply(babun.sens1, "[[", 2)) # effect estimate
estimate5.lb <- as.numeric(sapply(babun.sens1, "[[", 4)) # lower bound
estimate5.ub <- as.numeric(sapply(babun.sens1, "[[", 5)) # upper bound
pubs.out <- sapply(babun.sens1, "[[", 7) # publication

scores.list[[5]] <- (ma05$b - mean(estimate5))/sd(estimate5)

mean(estimate5.lb); mean(estimate5.ub)

# Richness: Carabids
ma06 <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=filter(brich, Taxon=="Carabids"), 
               method="REML",
               slab=EffectCode)

# jackknife analysis for abundance
brich.sens1 <- jackknife(filter(brich, Taxon=="Carabids"))
estimate6 <- as.numeric(sapply(brich.sens1, "[[", 2)) # effect estimate
estimate6.lb <- as.numeric(sapply(brich.sens1, "[[", 4)) # lower bound
estimate6.ub <- as.numeric(sapply(brich.sens1, "[[", 5)) # upper bound
pubs.out <- sapply(brich.sens1, "[[", 7) # publication

scores.list[[6]] <- (ma06$b - mean(estimate6))/sd(estimate6)

mean(estimate6.lb); mean(estimate6.ub)

# Bind all of the calculated effect sizes into a dataframe for plotting
TreatmentPlot <- as_tibble(rbind(c(8, "Abundance", "Bee Abundance", ma01$beta, ma01$ci.lb, ma01$ci.ub, ma01$pval, ma01$k),
                            c(7, "Richness", "Bee Richness", ma02$beta, ma02$ci.lb, ma02$ci.ub, ma02$pval, ma02$k),
                            c(2, "Abundance", "Butterfly Abundance", ma03$beta, ma03$ci.lb, ma03$ci.ub, ma03$pval, ma03$k),
                            #c(8, "Richness", "Butterflies", ma04$beta, ma04$ci.lb, ma04$ci.ub, ma04$pval, ma04$k),
                            c(5, "Abundance", "Carabidae Abundance", ma05$beta, ma05$ci.lb, ma05$ci.ub, ma05$pval, ma05$k),
                            c(4, "Richness", "Carabidae Richness", ma06$beta, ma06$ci.lb, ma06$ci.ub, ma06$pval, ma06$k)))

# Plot
ggplot(data=TreatmentPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
  geom_vline(xintercept=0, color="gray25", lty="dotted")+
  geom_errorbarh(aes(xmin=as.numeric(V5), xmax=as.numeric(V6)), height=0)+
  geom_point(size=3, pch=16, fill="black")+
  geom_text(aes(x=-4.0, label=paste(V3, " [",V8,"]", sep="")), hjust=0)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.line.x = element_line(color="black"))+
  xlab("Effect Size (g)")+
  ylim(1,9)+
  scale_x_continuous(limits=c(-5,5), breaks=c(-2, -1, -0.5, 0, 0.5, 1, 2))+
  ggtitle(label="Treatment+Fire Effects")

ggsave("TreatmentForestPlot.pdf", width=10, height=8, units="in")

# Specific Treatment Effects Richness
ma.c.gz <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=filter(babun, Taxon=="Carabids", Treatment=="Grazing"), 
               method="REML",
               slab=EffectCode)

ma.c.th.a <- rma.mv(yi=EffectSize, 
                V=EffectSizeStandardError,
                random=~1|PublicationCode/EffectCode,
                tdist=TRUE,
                data=filter(babun, Taxon=="Carabids", Treatment=="Thinning"), 
                method="REML",
                slab=EffectCode)

ma.c.th.r <- rma.mv(yi=EffectSize, 
                V=EffectSizeStandardError,
                random=~1|PublicationCode/EffectCode,
                tdist=TRUE,
                data=filter(brich, Taxon=="Carabids", Treatment=="Thinning"), 
                method="REML",
                slab=EffectCode)

ma.b.th.a <- rma.mv(yi=EffectSize, 
                V=EffectSizeStandardError,
                random=~1|PublicationCode/EffectCode,
                tdist=TRUE,
                data=filter(babun, Taxon=="Bees", Treatment=="Thinning"), 
                method="REML",
                slab=EffectCode)

ma.bu.th.a <- rma.mv(yi=EffectSize, 
                    V=EffectSizeStandardError,
                    random=~1|PublicationCode/EffectCode,
                    tdist=TRUE,
                    data=filter(babun, Taxon=="Butterflies", Treatment=="Thinning"), 
                    method="REML",
                    slab=EffectCode)

TreatmentPlot <- as_tibble(rbind(c(11, "Abundance", "Bee Abundance", ma01$beta, ma01$ci.lb, ma01$ci.ub, ma01$pval, ma01$k, "Major"),
                                 c(10, "Abundance", "Thinning + Fire", ma.b.th.a$beta, ma.b.th.a$ci.lb, ma.b.th.a$ci.ub, ma.b.th.a$pval, ma.b.th.a$k, "Minor"),
                                 c(12, "Richness", "Bee Richness", ma02$beta, ma02$ci.lb, ma02$ci.ub, ma02$pval, ma02$k, "Major"),
                                 c(2, "Abundance", "Butterfly Abundance", ma03$beta, ma03$ci.lb, ma03$ci.ub, ma03$pval, ma03$k, "Major"),
                                 c(1, "Abundance", "Thinning + Fire", ma.bu.th.a$beta, ma.bu.th.a$ci.lb, ma.bu.th.a$ci.ub, ma.bu.th.a$pval, ma.bu.th.a$k, "Minor"),
                                 #c(8, "Richness", "Butterflies", ma04$beta, ma04$ci.lb, ma04$ci.ub, ma04$pval, ma04$k),
                                 c(6, "Abundance", "Carabidae Abundance", ma05$beta, ma05$ci.lb, ma05$ci.ub, ma05$pval, ma05$k, "Major"),
                                 c(5, "Abundance", "Thinning + Fire", ma.c.th.a$beta, ma.c.th.a$ci.lb, ma.c.th.a$ci.ub, ma.c.th.a$pval, ma.c.th.a$k, "Minor"),
                                 c(4, "Abundance", "Grazing + Fire", ma.c.gz$beta, ma.c.gz$ci.lb, ma.c.gz$ci.ub, ma.c.gz$pval, ma.c.gz$k, "Minor"),
                                 c(8, "Richness", "Carabidae Richness", ma06$beta, ma06$ci.lb, ma06$ci.ub, ma06$pval, ma06$k, "Major"),
                                 c(7, "Richness", "Thinning + Fire", ma.c.th.r$beta, ma.c.th.r$ci.lb, ma.c.th.r$ci.ub, ma.c.th.r$pval, ma.c.th.r$k, "Minor")))

# Plot
ggplot(data=TreatmentPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
  geom_vline(xintercept=0, color="gray25", lty="dotted")+
  geom_errorbarh(aes(xmin=as.numeric(V5), xmax=as.numeric(V6)), height=0)+
  geom_point(size=3, pch=16, fill="black")+
  geom_text(data=dplyr::filter(TreatmentPlot, V9=="Minor"), aes(x=-3.0, label=paste(V3, " [",V8,"]", sep="")), hjust=0)+
  geom_text(data=dplyr::filter(TreatmentPlot, V9=="Major"), aes(x=-4.0, label=paste(V3, " [",V8,"]", sep="")), hjust=0)+
  geom_text(data=dplyr::filter(TreatmentPlot, V9=="Major"), aes(x=-4.0, label=paste(V3, " [",V8,"]", sep="")), hjust=0)+
  geom_text(data=dplyr::filter(TreatmentPlot, V9=="Major"), aes(x=-4.0, label=paste(V3, " [",V8,"]", sep="")), hjust=0)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.line.x = element_line(color="black"))+
  xlab("Effect Size (g)")+
  ylim(1,13)+
  scale_x_continuous(limits=c(-5,5), breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5))+
  ggtitle(label="Treatment+Fire Effects")

ggsave("TreatmentForestPlot_noColor.pdf", width=10, height=8, units="in")