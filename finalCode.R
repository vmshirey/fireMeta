## This script reads in effect size and standard error information from
## studies on insect responses to fire/other management practices,
## performs. 
##
##

# load required packages, establish re-utilized variables, custom functions
rm(list = ls())

require(tidyverse); require(bayesmeta); library(meta); library(ggmap); library(gridExtra); library(ggExtra); library(cowplot)
library(compute.es)
library(esc)

plots.list <- list() # list for storing plots
scores.list <- list() # list for storing z-scores

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

# load data containing papers, effect sizes + standard errors, and covariates
dat <- read.csv("FinalEffects.csv", stringsAsFactors = FALSE, sep=",", header=TRUE)
dat <- dat %>% filter(FireType != "Treatment+Fire", Author != "Warchola et al.")

# Filter by taxon of interest
dat <- dat %>% filter(Taxon=="Butterflies")

# subset data by taxon/diversity metric
babun <- dat %>% filter(Metric=="Abundance") # abundance
brich <- dat %>% filter(Metric=="Richness") # richness

###### RUN THE ANALYSIS ON THE FULL DATASET FOR ABUNDANCE METRIC #######
ma01 <- rma.mv(yi=EffectSize, 
               V=EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=babun, 
               method="REML",
               slab=EffectCode)

# jackknife analysis for abundance
babun.sens1 <- jackknife(babun)
estimate1 <- as.numeric(sapply(babun.sens1, "[[", 2)) # effect estimate
estimate1.lb <- as.numeric(sapply(babun.sens1, "[[", 4)) # lower bound
estimate1.ub <- as.numeric(sapply(babun.sens1, "[[", 5)) # upper bound
pubs.out <- sapply(babun.sens1, "[[", 7) # publication

scores.list[[1]] <- (ma01$b - mean(estimate1))/sd(estimate1)

# plot means, CIs for jackknife and actual effect
plots.list[[1]] <- qplot(x=seq(1,length(babun.sens1)), y=estimate1)+
    geom_errorbar(aes(x=seq(1,length(babun.sens1)), ymin=estimate1.lb, ymax=estimate1.ub, width=0.25))+
    geom_hline(yintercept=ma01$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=ma01$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=ma01$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate1), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate1.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate1.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens1)), y=estimate1))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon))

# forest plot for all studies and abundance
plots.list[[2]] <- forest(ma01, order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)

###### RUN THE ANALYSIS ON THE FULL DATASET FOR RICHNESS METRIC #######
ma02 <- rma.mv(yi=EffectSize, 
               V=brich$EffectSizeStandardError,
               random=~1|PublicationCode/EffectCode,
               tdist=TRUE,
               data=brich, 
               method="REML",
               slab=EffectCode)

# jackknife analysis for abundance
brich.sens1 <- jackknife(brich)
estimate2 <- as.numeric(sapply(brich.sens1, "[[", 2))
estimate2.lb <- as.numeric(sapply(brich.sens1, "[[", 4))
estimate2.ub <- as.numeric(sapply(brich.sens1, "[[", 5))
pubs.out <- sapply(brich.sens1, "[[", 7)

scores.list[[2]] <- (ma02$b - mean(estimate2))/sd(estimate2)

# plot mean variation
boxplot(estimate2)
abline(h=ma02$b, col="#D55E00", lwd=2, lty=2)

# plot means, CIs for jackknife and actual effect
plots.list[[3]] <- qplot(x=seq(1,length(brich.sens1)), y=estimate2)+
    geom_errorbar(aes(x=seq(1,length(brich.sens1)), ymin=estimate2.lb, ymax=estimate2.ub, width=0.25))+
    geom_hline(yintercept=ma02$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=ma02$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=ma02$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate2), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate2.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate2.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens1)), y=estimate2))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon))

# forest plot for all studies and richness
plots.list[[4]] <- forest(ma02, order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)

# Fire Severity on Abundance
babun.hi <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode, subset=(Severity=="High"), slab=EffectCode)
babun.mi <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode, subset=(Severity=="Mid"), slab=EffectCode)
babun.lo <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode, subset=(Severity=="Low"), slab=EffectCode)

babun.sens.hi <- jackknife(subset(babun, Severity=="High"))
estimate3 <- as.numeric(sapply(babun.sens.hi, "[[", 2))
estimate3.lb <- as.numeric(sapply(babun.sens.hi, "[[", 4))
estimate3.ub <- as.numeric(sapply(babun.sens.hi, "[[", 5))
pubs.out <- sapply(babun.sens.hi, "[[", 7)

scores.list[[3]] <- (babun.hi$b - mean(estimate3))/sd(estimate3)

plots.list[[4]] <- qplot(x=seq(1,length(babun.sens.hi)), y=estimate3)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.hi)), ymin=estimate3.lb, ymax=estimate3.ub, width=0.25))+
    geom_hline(yintercept=babun.hi$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.hi$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.hi$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate3), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate3.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate3.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.hi)), y=estimate3))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "High Severity"))

babun.sens.mi <- jackknife(subset(babun, Severity=="Mid"))
estimate4 <- as.numeric(sapply(babun.sens.mi, "[[", 2))
estimate4.lb <- as.numeric(sapply(babun.sens.mi, "[[", 4))
estimate4.ub <- as.numeric(sapply(babun.sens.mi, "[[", 5))
pubs.out <- sapply(babun.sens.mi, "[[", 7)

scores.list[[4]] <- (babun.mi$b - mean(estimate4))/sd(estimate4)

plots.list[[5]] <- qplot(x=seq(1,length(babun.sens.mi)), y=estimate4)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.mi)), ymin=estimate4.lb, ymax=estimate4.ub, width=0.25))+
    geom_hline(yintercept=babun.mi$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.mi$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.mi$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate4), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate4.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate4.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.mi)), y=estimate4))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Mid Severity"))

babun.sens.lo <- jackknife(subset(babun, Severity=="Low"))
estimate5 <- as.numeric(sapply(babun.sens.lo, "[[", 2))
estimate5.lb <- as.numeric(sapply(babun.sens.lo, "[[", 4))
estimate5.ub <- as.numeric(sapply(babun.sens.lo, "[[", 5))
pubs.out <- sapply(babun.sens.lo, "[[", 7)

scores.list[[5]] <- (babun.lo$b - mean(estimate5))/sd(estimate5)

plots.list[[6]] <-qplot(x=seq(1,length(babun.sens.lo)), y=estimate5)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.lo)), ymin=estimate5.lb, ymax=estimate5.ub, width=0.25))+
    geom_hline(yintercept=babun.lo$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.lo$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.lo$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate5), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate5.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate5.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.lo)), y=estimate5))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Low Severity"))


# Fire Severity on Richness
brich.hi <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(Severity=="High"), slab=EffectCode)
brich.mi <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(Severity=="Mid"), slab=EffectCode)
brich.lo <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(Severity=="Low"), slab=EffectCode)

brich.sens.hi <- jackknife(subset(brich, Severity=="High"))
estimate6 <- as.numeric(sapply(brich.sens.hi, "[[", 2))
estimate6.lb <- as.numeric(sapply(brich.sens.hi, "[[", 4))
estimate6.ub <- as.numeric(sapply(brich.sens.hi, "[[", 5))
pubs.out <- sapply(brich.sens.hi, "[[", 7)

scores.list[[6]] <- (brich.hi$b - mean(estimate6))/sd(estimate6)

plots.list[[7]] <- qplot(x=seq(1,length(brich.sens.hi)), y=estimate6)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.hi)), ymin=estimate6.lb, ymax=estimate6.ub, width=0.25))+
    geom_hline(yintercept=brich.hi$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.hi$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.hi$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate6), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate6.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate6.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.hi)), y=estimate6))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "High Severity"))

brich.sens.mi <- jackknife(subset(brich, Severity=="Mid"))
estimate7 <- as.numeric(sapply(brich.sens.mi, "[[", 2))
estimate7.lb <- as.numeric(sapply(brich.sens.mi, "[[", 4))
estimate7.ub <- as.numeric(sapply(brich.sens.mi, "[[", 5))
pubs.out <- sapply(brich.sens.mi, "[[", 7)

scores.list[[7]] <- (brich.mi$b - mean(estimate7))/sd(estimate7)

plots.list[[8]] <-qplot(x=seq(1,length(brich.sens.mi)), y=estimate7)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.mi)), ymin=estimate7.lb, ymax=estimate7.ub, width=0.25))+
    geom_hline(yintercept=brich.mi$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.mi$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.mi$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate7), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate7.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate7.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.mi)), y=estimate7))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Mid Severity"))

brich.sens.lo <- jackknife(subset(brich, Severity=="Low"))
estimate8 <- as.numeric(sapply(brich.sens.lo, "[[", 2))
estimate8.lb <- as.numeric(sapply(brich.sens.lo, "[[", 4))
estimate8.ub <- as.numeric(sapply(brich.sens.lo, "[[", 5))
pubs.out <- sapply(brich.sens.lo, "[[", 7)

scores.list[[8]] <- (brich.lo$b - mean(estimate8))/sd(estimate8)

plots.list[[9]] <-qplot(x=seq(1,length(brich.sens.lo)), y=estimate8)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.lo)), ymin=estimate8.lb, ymax=estimate8.ub, width=0.25))+
    geom_hline(yintercept=brich.lo$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.lo$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.lo$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate8), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate8.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate8.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.lo)), y=estimate8))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Low Severity"))

# Growing vs. Non-growing On Abundance
babun.g <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(Growing.non=="Growing"), slab=EffectCode)
babun.ng <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(Growing.non=="Dormant"), slab=EffectCode)
babun.trp <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(Growing.non=="Tropics"), slab=EffectCode)

babun.sens.g <- jackknife(subset(babun, Growing.non=="Growing"))
estimate9 <- as.numeric(sapply(babun.sens.g, "[[", 2))
estimate9.lb <- as.numeric(sapply(babun.sens.g, "[[", 4))
estimate9.ub <- as.numeric(sapply(babun.sens.g, "[[", 5))
pubs.out <- sapply(babun.sens.g, "[[", 7)

scores.list[[9]] <- (babun.g$b - mean(estimate9))/sd(estimate9)

plots.list[[10]] <- qplot(x=seq(1,length(babun.sens.g)), y=estimate9)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.g)), ymin=estimate9.lb, ymax=estimate9.ub, width=0.25))+
    geom_hline(yintercept=babun.g$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.g$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.g$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate9), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate9.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate9.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.g)), y=estimate9))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Growing Season"))

babun.sens.ng <- jackknife(subset(babun, Growing.non=="Dormant"))
estimate10 <- as.numeric(sapply(babun.sens.ng, "[[", 2))
estimate10.lb <- as.numeric(sapply(babun.sens.ng, "[[", 4))
estimate10.ub <- as.numeric(sapply(babun.sens.ng, "[[", 5))
pubs.out <- sapply(babun.sens.ng, "[[", 7)

scores.list[[10]] <- (babun.ng$b - mean(estimate10))/sd(estimate10)

plots.list[[11]] <- qplot(x=seq(1,length(babun.sens.ng)), y=estimate10)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.ng)), ymin=estimate10.lb, ymax=estimate10.ub, width=0.25))+
    geom_hline(yintercept=babun.ng$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.ng$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.ng$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate10), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate10.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate10.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.ng)), y=estimate10))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Non-Growing Season"))

# Growing vs. Non-growing On Richness
brich.g <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(Growing.non=="Growing"), slab=EffectCode)
brich.ng <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(Growing.non=="Dormant"), slab=EffectCode)

brich.sens.g <- jackknife(subset(brich, Growing.non=="Growing"))
estimate11 <- as.numeric(sapply(brich.sens.g, "[[", 2))
estimate11.lb <- as.numeric(sapply(brich.sens.g, "[[", 4))
estimate11.ub <- as.numeric(sapply(brich.sens.g, "[[", 5))
pubs.out <- sapply(brich.sens.g, "[[", 7)

scores.list[[10]] <- (brich.g$b - mean(estimate11))/sd(estimate11)

plots.list[[12]] <- qplot(x=seq(1,length(brich.sens.g)), y=estimate11)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.g)), ymin=estimate11.lb, ymax=estimate11.ub, width=0.25))+
    geom_hline(yintercept=brich.g$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.g$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.g$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate11), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate11.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate11.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.g)), y=estimate11))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Growing Season"))

brich.sens.ng <- jackknife(subset(brich, Growing.non=="Dormant"))
estimate12 <- as.numeric(sapply(brich.sens.ng, "[[", 2))
estimate12.lb <- as.numeric(sapply(brich.sens.ng, "[[", 4))
estimate12.ub <- as.numeric(sapply(brich.sens.ng, "[[", 5))
pubs.out <- sapply(brich.sens.ng, "[[", 7)

scores.list[[11]] <- (brich.ng$b - mean(estimate12))/sd(estimate12)

plots.list[[13]] <- qplot(x=seq(1,length(brich.sens.ng)), y=estimate12)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.ng)), ymin=estimate12.lb, ymax=estimate12.ub, width=0.25))+
    geom_hline(yintercept=brich.ng$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.ng$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.ng$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate12), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate12.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate12.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.ng)), y=estimate12))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Non-Growing Season"))

# Prescribed vs. Wildfire vs. Treatment+Fire on Abundance
babun.pf <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireType=="Prescribed Fire"), slab=EffectCode)
babun.wf <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireType=="Wildfire"), slab=EffectCode)
# babun.tf <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireType=="Treatment+Fire"), slab=EffectCode)

babun.sens.pf <- jackknife(subset(babun, FireType=="Prescribed Fire"))
estimate13 <- as.numeric(sapply(babun.sens.pf, "[[", 2))
estimate13.lb <- as.numeric(sapply(babun.sens.pf, "[[", 4))
estimate13.ub <- as.numeric(sapply(babun.sens.pf, "[[", 5))
pubs.out <- sapply(babun.sens.pf, "[[", 7)

scores.list[[12]] <- (babun.pf$b - mean(estimate13))/sd(estimate13)

plots.list[[14]] <- qplot(x=seq(1,length(babun.sens.pf)), y=estimate13)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.pf)), ymin=estimate13.lb, ymax=estimate13.ub, width=0.25))+
    geom_hline(yintercept=babun.pf$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.pf$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.pf$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate13), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate13.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate13.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.pf)), y=estimate13))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Prescribed Fire"))

babun.sens.wf <- jackknife(subset(babun, FireType=="Wildfire"))
estimate14 <- as.numeric(sapply(babun.sens.wf, "[[", 2))
estimate14.lb <- as.numeric(sapply(babun.sens.wf, "[[", 4))
estimate14.ub <- as.numeric(sapply(babun.sens.wf, "[[", 5))
pubs.out <- sapply(babun.sens.wf, "[[", 7)

scores.list[[13]] <- (babun.wf$b - mean(estimate14))/sd(estimate14)

plots.list[[15]] <- qplot(x=seq(1,length(babun.sens.wf)), y=estimate14)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.wf)), ymin=estimate14.lb, ymax=estimate14.ub, width=0.25))+
    geom_hline(yintercept=babun.wf$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.wf$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.wf$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate14), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate14.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate14.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.wf)), y=estimate14))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Wildfire"))

# babun.sens.tf <- jackknife(subset(babun, FireType=="Treatment+Fire"))
# estimate15 <- as.numeric(sapply(babun.sens.tf, "[[", 2))
# estimate15.lb <- as.numeric(sapply(babun.sens.tf, "[[", 4))
# estimate15.ub <- as.numeric(sapply(babun.sens.tf, "[[", 5))
# pubs.out <- sapply(babun.sens.tf, "[[", 7)
# 
# scores.list[[14]] <- (babun.tf$b - mean(estimate15))/sd(estimate15)
# 
# plots.list[[16]] <- qplot(x=seq(1,length(babun.sens.tf)), y=estimate15)+
#     geom_errorbar(aes(x=seq(1,length(babun.sens.tf)), ymin=estimate15.lb, ymax=estimate15.ub, width=0.25))+
#     geom_hline(yintercept=babun.tf$b, color="#D55E00", size=1.5)+
#     geom_hline(yintercept=babun.tf$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
#     geom_hline(yintercept=babun.tf$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
#     geom_hline(yintercept=mean(estimate15), color="#0072B2", size=1.5)+
#     geom_hline(yintercept=mean(estimate15.lb), color="#0072B2", linetype="dashed", size=1.5)+
#     geom_hline(yintercept=mean(estimate15.ub), color="#0072B2", linetype="dashed", size=1.5)+
#     geom_point(aes(x=seq(1,length(babun.sens.tf)), y=estimate15))+
#     ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Treatment + Fire"))

# Prescribed vs. Wildfire vs. Treatment+Fire on Richness
brich.pf <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireType=="Prescribed Fire"), slab=EffectCode)
brich.wf <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireType=="Wildfire"), slab=EffectCode)
# brich.tf <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireType=="Treatment+Fire"), slab=EffectCode)

brich.sens.pf <- jackknife(subset(brich, FireType=="Prescribed Fire"))
estimate16 <- as.numeric(sapply(brich.sens.pf, "[[", 2))
estimate16.lb <- as.numeric(sapply(brich.sens.pf, "[[", 4))
estimate16.ub <- as.numeric(sapply(brich.sens.pf, "[[", 5))
pubs.out <- sapply(brich.sens.pf, "[[", 7)

scores.list[[15]] <- (brich.pf$b - mean(estimate16))/sd(estimate16)

plots.list[[17]] <- qplot(x=seq(1,length(brich.sens.pf)), y=estimate16)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.pf)), ymin=estimate16.lb, ymax=estimate16.ub, width=0.25))+
    geom_hline(yintercept=brich.pf$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.pf$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.pf$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate16), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate16.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate16.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.pf)), y=estimate16))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Prescribed Fire"))

brich.sens.wf <- jackknife(subset(brich, FireType=="Wildfire"))
estimate17 <- as.numeric(sapply(brich.sens.wf, "[[", 2))
estimate17.lb <- as.numeric(sapply(brich.sens.wf, "[[", 4))
estimate17.ub <- as.numeric(sapply(brich.sens.wf, "[[", 5))
pubs.out <- sapply(brich.sens.wf, "[[", 7)

scores.list[[16]] <- (brich.wf$b - mean(estimate17))/sd(estimate17)

plots.list[[18]] <- qplot(x=seq(1,length(brich.sens.wf)), y=estimate17)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.wf)), ymin=estimate17.lb, ymax=estimate17.ub, width=0.25))+
    geom_hline(yintercept=brich.wf$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.wf$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.wf$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate17), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate17.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate17.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.wf)), y=estimate17))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Wildfire"))

# brich.sens.tf <- jackknife(subset(brich, FireType=="Treatment+Fire"))
# estimate18 <- as.numeric(sapply(brich.sens.tf, "[[", 2))
# estimate18.lb <- as.numeric(sapply(brich.sens.tf, "[[", 4))
# estimate18.ub <- as.numeric(sapply(brich.sens.tf, "[[", 5))
# pubs.out <- sapply(brich.sens.tf, "[[", 7)
# 
# scores.list[[17]] <- (brich.tf$b - mean(estimate18))/sd(estimate18)

# plots.list[[19]] <- qplot(x=seq(1,length(brich.sens.tf)), y=estimate18)+
#     geom_errorbar(aes(x=seq(1,length(brich.sens.tf)), ymin=estimate18.lb, ymax=estimate18.ub, width=0.25))+
#     geom_hline(yintercept=brich.tf$b, color="#D55E00", size=1.5)+
#     geom_hline(yintercept=brich.tf$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
#     geom_hline(yintercept=brich.tf$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
#     geom_hline(yintercept=mean(estimate18), color="#0072B2", size=1.5)+
#     geom_hline(yintercept=mean(estimate18.lb), color="#0072B2", linetype="dashed", size=1.5)+
#     geom_hline(yintercept=mean(estimate18.ub), color="#0072B2", linetype="dashed", size=1.5)+
#     geom_point(aes(x=seq(1,length(brich.sens.tf)), y=estimate18))+
#     ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Treatment + Fire"))

# Time since fire (0, 1to5, 5to10, 10+) on Abundance
babun.ze <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(TimeClass=="0"), slab=EffectCode)
babun.of <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(TimeClass=="1to5"), slab=EffectCode)
babun.ft <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(TimeClass=="5to10"), slab=EffectCode)
babun.tp <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(TimeClass=="10+"), slab=EffectCode)

babun.sens.ze <- jackknife(subset(babun, TimeClass=="0"))
estimate19 <- as.numeric(sapply(babun.sens.ze, "[[", 2))
estimate19.lb <- as.numeric(sapply(babun.sens.ze, "[[", 4))
estimate19.ub <- as.numeric(sapply(babun.sens.ze, "[[", 5))
pubs.out <- sapply(babun.sens.ze, "[[", 7)

scores.list[[18]] <- (babun.ze$b - mean(estimate19))/sd(estimate19)

plots.list[[20]] <- qplot(x=seq(1,length(babun.sens.ze)), y=estimate19)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.ze)), ymin=estimate19.lb, ymax=estimate19.ub, width=0.25))+
    geom_hline(yintercept=babun.ze$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.ze$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.ze$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate19), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate19.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate19.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.ze)), y=estimate19))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Zero years"))

babun.sens.of <- jackknife(subset(babun, TimeClass=="1to5"))
estimate20 <- as.numeric(sapply(babun.sens.of, "[[", 2))
estimate20.lb <- as.numeric(sapply(babun.sens.of, "[[", 4))
estimate20.ub <- as.numeric(sapply(babun.sens.of, "[[", 5))
pubs.out <- sapply(babun.sens.of, "[[", 7)

scores.list[[19]] <- (babun.of$b - mean(estimate20))/sd(estimate20)

plots.list[[21]] <- qplot(x=seq(1,length(babun.sens.of)), y=estimate20)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.of)), ymin=estimate20.lb, ymax=estimate20.ub, width=0.25))+
    geom_hline(yintercept=babun.of$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.of$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.of$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate20), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate20.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate20.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.of)), y=estimate20))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "1-5 years"))

babun.sens.ft <- jackknife(subset(babun, TimeClass=="5to10"))
estimate21 <- as.numeric(sapply(babun.sens.ft, "[[", 2))
estimate21.lb <- as.numeric(sapply(babun.sens.ft, "[[", 4))
estimate21.ub <- as.numeric(sapply(babun.sens.ft, "[[", 5))
pubs.out <- sapply(babun.sens.ft, "[[", 7)

scores.list[[20]] <- (babun.ft$b - mean(estimate21))/sd(estimate21)

plots.list[[22]] <- qplot(x=seq(1,length(babun.sens.ft)), y=estimate21)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.ft)), ymin=estimate21.lb, ymax=estimate21.ub, width=0.25))+
    geom_hline(yintercept=babun.ft$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.ft$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.ft$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate21), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate21.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate21.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.ft)), y=estimate21))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "5-10 years"))

babun.sens.tp <- jackknife(subset(babun, TimeClass=="10+"))
estimate22 <- as.numeric(sapply(babun.sens.tp, "[[", 2))
estimate22.lb <- as.numeric(sapply(babun.sens.tp, "[[", 4))
estimate22.ub <- as.numeric(sapply(babun.sens.tp, "[[", 5))
pubs.out <- sapply(babun.sens.tp, "[[", 7)

scores.list[[21]] <- (babun.tp$b - mean(estimate22))/sd(estimate22)

plots.list[[23]] <- qplot(x=seq(1,length(babun.sens.tp)), y=estimate22)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.tp)), ymin=estimate22.lb, ymax=estimate22.ub, width=0.25))+
    geom_hline(yintercept=babun.tp$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.tp$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.tp$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate22), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate22.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate22.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.tp)), y=estimate22))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "10+ years"))

# Time since fire (0, 1to5, 5to10, 10+) on Richness
brich.ze <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(TimeClass=="0"), slab=EffectCode)
brich.of <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(TimeClass=="1to5"), slab=EffectCode)
brich.ft <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(TimeClass=="5to10"), slab=EffectCode)
brich.tp <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(TimeClass=="10+"), slab=EffectCode)

brich.sens.ze <- jackknife(subset(brich, TimeClass=="0"))
estimate23 <- as.numeric(sapply(brich.sens.ze, "[[", 2))
estimate23.lb <- as.numeric(sapply(brich.sens.ze, "[[", 4))
estimate23.ub <- as.numeric(sapply(brich.sens.ze, "[[", 5))
pubs.out <- sapply(brich.sens.ze, "[[", 7)

scores.list[[22]] <- (brich.ze$b - mean(estimate23))/sd(estimate23)

plots.list[[24]] <- qplot(x=seq(1,length(brich.sens.ze)), y=estimate23)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.ze)), ymin=estimate23.lb, ymax=estimate23.ub, width=0.25))+
    geom_hline(yintercept=brich.ze$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.ze$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.ze$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate23), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate23.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate23.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.ze)), y=estimate23))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Zero years"))

brich.sens.of <- jackknife(subset(brich, TimeClass=="1to5"))
estimate24 <- as.numeric(sapply(brich.sens.of, "[[", 2))
estimate24.lb <- as.numeric(sapply(brich.sens.of, "[[", 4))
estimate24.ub <- as.numeric(sapply(brich.sens.of, "[[", 5))
pubs.out <- sapply(brich.sens.of, "[[", 7)

scores.list[[23]] <- (brich.of$b - mean(estimate24))/sd(estimate24)

plots.list[[25]] <- qplot(x=seq(1,length(brich.sens.of)), y=estimate24)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.of)), ymin=estimate24.lb, ymax=estimate24.ub, width=0.25))+
    geom_hline(yintercept=brich.of$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.of$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.of$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate24), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate24.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate24.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.of)), y=estimate24))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "1-5 years"))

brich.sens.ft <- jackknife(subset(brich, TimeClass=="5to10"))
estimate25 <- as.numeric(sapply(brich.sens.ft, "[[", 2))
estimate25.lb <- as.numeric(sapply(brich.sens.ft, "[[", 4))
estimate25.ub <- as.numeric(sapply(brich.sens.ft, "[[", 5))
pubs.out <- sapply(brich.sens.ft, "[[", 7)

scores.list[[24]] <- (brich.ft$b - mean(estimate25))/sd(estimate25)

plots.list[[26]] <- qplot(x=seq(1,length(brich.sens.ft)), y=estimate25)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.ft)), ymin=estimate25.lb, ymax=estimate25.ub, width=0.25))+
    geom_hline(yintercept=brich.ft$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.ft$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.ft$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate25), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate25.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate25.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.ft)), y=estimate25))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "5-10 years"))

brich.sens.tp <- jackknife(subset(brich, TimeClass=="10+"))
estimate26 <- as.numeric(sapply(brich.sens.tp, "[[", 2))
estimate26.lb <- as.numeric(sapply(brich.sens.tp, "[[", 4))
estimate26.ub <- as.numeric(sapply(brich.sens.tp, "[[", 5))
pubs.out <- sapply(brich.sens.tp, "[[", 7)

scores.list[[25]] <- (brich.tp$b - mean(estimate26))/sd(estimate26)

plots.list[[27]] <- qplot(x=seq(1,length(brich.sens.tp)), y=estimate26)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.tp)), ymin=estimate26.lb, ymax=estimate26.ub, width=0.25))+
    geom_hline(yintercept=brich.tp$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.tp$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.tp$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate26), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate26.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate26.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.tp)), y=estimate26))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "10+ years"))

# Fire Freq (yes/no) On Abundance
babun.ff <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireFreq=="Yes"), slab=EffectCode)
babun.nff <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireFreq=="No"), slab=EffectCode)

babun.sens.ff <- jackknife(subset(babun, FireFreq=="Yes"))
estimate27 <- as.numeric(sapply(babun.sens.ff, "[[", 2))
estimate27.lb <- as.numeric(sapply(babun.sens.ff, "[[", 4))
estimate27.ub <- as.numeric(sapply(babun.sens.ff, "[[", 5))
pubs.out <- sapply(babun.sens.ff, "[[", 7)

scores.list[[26]] <- (babun.ff$b - mean(estimate27))/sd(estimate27)

plots.list[[28]] <- qplot(x=seq(1,length(babun.sens.ff)), y=estimate27)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.ff)), ymin=estimate27.lb, ymax=estimate27.ub, width=0.25))+
    geom_hline(yintercept=babun.ff$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.ff$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.ff$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate27), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate27.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate27.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.ff)), y=estimate27))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Fire Frequency"))

babun.sens.nff <- jackknife(subset(babun, FireFreq=="No"))
estimate28 <- as.numeric(sapply(babun.sens.nff, "[[", 2))
estimate28.lb <- as.numeric(sapply(babun.sens.nff, "[[", 4))
estimate28.ub <- as.numeric(sapply(babun.sens.nff, "[[", 5))
pubs.out <- sapply(babun.sens.nff, "[[", 7)

scores.list[[27]] <- (babun.nff$b - mean(estimate28))/sd(estimate28)

plots.list[[29]] <- qplot(x=seq(1,length(babun.sens.nff)), y=estimate28)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.nff)), ymin=estimate28.lb, ymax=estimate28.ub, width=0.25))+
    geom_hline(yintercept=babun.nff$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.nff$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.nff$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate28), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate28.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate28.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.nff)), y=estimate28))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "No Fire Frequency"))

# Fire Freq (yes/no) On Richness
brich.ff <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireFreq=="Yes"), slab=EffectCode)
brich.nff <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(FireFreq=="No"), slab=EffectCode)

brich.sens.ff <- jackknife(subset(brich, Growing.non=="Yes"))
estimate29 <- as.numeric(sapply(brich.sens.ff, "[[", 2))
estimate29.lb <- as.numeric(sapply(brich.sens.ff, "[[", 4))
estimate29.ub <- as.numeric(sapply(brich.sens.ff, "[[", 5))
pubs.out <- sapply(brich.sens.ff, "[[", 7)

scores.list[[28]] <- (brich.ff$b - mean(estimate29))/sd(estimate29)

plots.list[[30]] <- qplot(x=seq(1,length(brich.sens.ff)), y=estimate29)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.ff)), ymin=estimate29.lb, ymax=estimate29.ub, width=0.25))+
    geom_hline(yintercept=brich.ff$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.ff$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.ff$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate29), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate29.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate29.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.ff)), y=estimate29))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Fire Frequency"))

brich.sens.nff <- jackknife(subset(brich, FireFreq=="No"))
estimate30 <- as.numeric(sapply(brich.sens.nff, "[[", 2))
estimate30.lb <- as.numeric(sapply(brich.sens.nff, "[[", 4))
estimate30.ub <- as.numeric(sapply(brich.sens.nff, "[[", 5))
pubs.out <- sapply(brich.sens.nff, "[[", 7)

scores.list[[29]] <- (brich.nff$b - mean(estimate30))/sd(estimate30)

plots.list[[31]] <- qplot(x=seq(1,length(brich.sens.nff)), y=estimate30)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.nff)), ymin=estimate30.lb, ymax=estimate30.ub, width=0.25))+
    geom_hline(yintercept=brich.nff$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.nff$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.nff$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate30), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate30.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate30.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.nff)), y=estimate30))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "No Fire Frequency"))

# Habitat (Deciduous Forest, Grassland, Coniferous Forest) on Abundance
babun.df <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=list(~1|PublicationCode,                             
                                                                                                    ~1|EffectCode),subset=(HabitatType=="Broad-Leaf Forest"), slab=EffectCode)
babun.cf <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=list(~1|PublicationCode,                             
                                                                                                   ~1|EffectCode),subset=(HabitatType=="Coniferous Forest"), slab=EffectCode)
babun.gr <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", random=list(~1|PublicationCode,                              
                                                                                                   ~1|EffectCode),subset=(HabitatType=="Grassland"), slab=EffectCode)

babun.sens.df <- jackknife(subset(babun, HabitatType=="Broad-Leaf Forest"))
estimate31 <- as.numeric(sapply(babun.sens.df, "[[", 2))
estimate31.lb <- as.numeric(sapply(babun.sens.df, "[[", 4))
estimate31.ub <- as.numeric(sapply(babun.sens.df, "[[", 5))
pubs.out <- sapply(babun.sens.df, "[[", 7)

scores.list[[30]] <- (babun.df$b - mean(estimate31))/sd(estimate31)

plots.list[[32]] <-  qplot(x=seq(1,length(babun.sens.df)), y=estimate31)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.df)), ymin=estimate31.lb, ymax=estimate31.ub, width=0.25))+
    geom_hline(yintercept=babun.df$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.df$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.df$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate31), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate31.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate31.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.df)), y=estimate31))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Broad-Leaf Forest"))

babun.sens.cf <- jackknife(subset(babun, HabitatType=="Coniferous Forest"))
estimate32 <- as.numeric(sapply(babun.sens.cf, "[[", 2))
estimate32.lb <- as.numeric(sapply(babun.sens.cf, "[[", 4))
estimate32.ub <- as.numeric(sapply(babun.sens.cf, "[[", 5))
pubs.out <- sapply(babun.sens.cf, "[[", 7)

scores.list[[31]] <- (babun.cf$b - mean(estimate32))/sd(estimate32)

plots.list[[33]] <-  qplot(x=seq(1,length(babun.sens.cf)), y=estimate32)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.cf)), ymin=estimate32.lb, ymax=estimate32.ub, width=0.25))+
    geom_hline(yintercept=babun.cf$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.cf$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.cf$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate32), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate32.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate32.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.cf)), y=estimate32))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Coniferous Forest"))

babun.sens.gr <- jackknife(subset(babun, HabitatType=="Grassland"))
estimate33 <- as.numeric(sapply(babun.sens.gr, "[[", 2))
estimate33.lb <- as.numeric(sapply(babun.sens.gr, "[[", 4))
estimate33.ub <- as.numeric(sapply(babun.sens.gr, "[[", 5))
pubs.out <- sapply(babun.sens.gr, "[[", 7)

scores.list[[32]] <- (babun.gr$b - mean(estimate33))/sd(estimate33)

plots.list[[34]] <-  qplot(x=seq(1,length(babun.sens.gr)), y=estimate33)+
    geom_errorbar(aes(x=seq(1,length(babun.sens.gr)), ymin=estimate33.lb, ymax=estimate33.ub, width=0.25))+
    geom_hline(yintercept=babun.gr$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=babun.gr$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=babun.gr$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate33), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate33.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate33.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(babun.sens.gr)), y=estimate33))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon, "Grasslands"))

# Habitat (Deciduous Forest, Grassland, Coniferous Forest) on Richness
brich.df <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(HabitatType=="Broad-Leaf Forest"), slab=EffectCode)
brich.cf <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(HabitatType=="Coniferous Forest"), slab=EffectCode)
brich.gr <-rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", random=~1|PublicationCode/EffectCode,subset=(HabitatType=="Grassland"), slab=EffectCode)

brich.sens.df <- jackknife(subset(brich, HabitatType=="Broad-Leaf Forest"))
estimate34 <- as.numeric(sapply(brich.sens.df, "[[", 2))
estimate34.lb <- as.numeric(sapply(brich.sens.df, "[[", 4))
estimate34.ub <- as.numeric(sapply(brich.sens.df, "[[", 5))
pubs.out <- sapply(brich.sens.df, "[[", 7)

scores.list[[33]] <- (brich.df$b - mean(estimate34))/sd(estimate34)

plots.list[[35]] <- qplot(x=seq(1,length(brich.sens.df)), y=estimate34)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.df)), ymin=estimate34.lb, ymax=estimate34.ub, width=0.25))+
    geom_hline(yintercept=brich.df$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.df$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.df$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate34), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate34.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate34.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.df)), y=estimate34))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Broad-Leaf Forest"))

brich.sens.cf <- jackknife(subset(brich, HabitatType=="Coniferous Forest"))
estimate35 <- as.numeric(sapply(brich.sens.cf, "[[", 2))
estimate35.lb <- as.numeric(sapply(brich.sens.cf, "[[", 4))
estimate35.ub <- as.numeric(sapply(brich.sens.cf, "[[", 5))
pubs.out <- sapply(brich.sens.cf, "[[", 7)

scores.list[[34]] <- (brich.cf$b - mean(estimate35))/sd(estimate35)

plots.list[[36]] <- qplot(x=seq(1,length(brich.sens.cf)), y=estimate35)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.cf)), ymin=estimate35.lb, ymax=estimate35.ub, width=0.25))+
    geom_hline(yintercept=brich.cf$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.cf$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.cf$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate35), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate35.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate35.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.cf)), y=estimate35))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Coniferous Forest"))

brich.sens.gr <- jackknife(subset(brich, HabitatType=="Grassland"))
estimate36 <- as.numeric(sapply(brich.sens.gr, "[[", 2))
estimate36.lb <- as.numeric(sapply(brich.sens.gr, "[[", 4))
estimate36.ub <- as.numeric(sapply(brich.sens.gr, "[[", 5))
pubs.out <- sapply(brich.sens.gr, "[[", 7)

scores.list[[35]] <- (brich.gr$b - mean(estimate36))/sd(estimate36)

plots.list[[37]] <- qplot(x=seq(1,length(brich.sens.gr)), y=estimate36)+
    geom_errorbar(aes(x=seq(1,length(brich.sens.gr)), ymin=estimate36.lb, ymax=estimate36.ub, width=0.25))+
    geom_hline(yintercept=brich.gr$b, color="#D55E00", size=1.5)+
    geom_hline(yintercept=brich.gr$ci.lb, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=brich.gr$ci.ub, color="#D55E00", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate36), color="#0072B2", size=1.5)+
    geom_hline(yintercept=mean(estimate36.lb), color="#0072B2", linetype="dashed", size=1.5)+
    geom_hline(yintercept=mean(estimate36.ub), color="#0072B2", linetype="dashed", size=1.5)+
    geom_point(aes(x=seq(1,length(brich.sens.gr)), y=estimate36))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon, "Grasslands"))

# Funnel plots
funnel.plots <- list()

funnel.plots[[1]] <- funnel(ma01, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[2]] <- funnel(ma02, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[3]] <- funnel(babun.wf, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[4]] <- funnel(babun.pf, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
# funnel.plots[[5]] <- funnel(babun.tf, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[6]] <- funnel(brich.wf, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[7]] <- funnel(brich.pf, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
# funnel.plots[[8]] <- funnel(brich.tf, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[9]] <- funnel(babun.g, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[10]] <- funnel(babun.ng, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[11]] <- funnel(brich.g, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[12]] <- funnel(brich.ng, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[13]] <- funnel(babun.df, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[14]] <- funnel(babun.cf, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[15]] <- funnel(babun.gr, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[16]] <- funnel(brich.df, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[17]] <- funnel(brich.cf, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[18]] <- funnel(brich.gr, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[19]] <- funnel(babun.lo, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[20]] <- funnel(babun.mi, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[21]] <- funnel(babun.hi, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[22]] <- funnel(brich.lo, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[23]] <- funnel(brich.mi, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[24]] <- funnel(brich.hi, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[25]] <- funnel(babun.ze, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[26]] <- funnel(babun.of, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[27]] <- funnel(babun.ft, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[28]] <- funnel(babun.tp, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[29]] <- funnel(brich.ze, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[30]] <- funnel(brich.of, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[31]] <- funnel(brich.ft, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[32]] <- funnel(brich.tp, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[33]] <- funnel(babun.ff, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[34]] <- funnel(babun.nff, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

funnel.plots[[35]] <- funnel(brich.ff, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))
funnel.plots[[36]] <- funnel(brich.nff, legend=TRUE, level=c(90, 95, 99), shade=c("white", "#EE442F", "#63ACBE"))

#### FOREST PLOT ####
# Abundance
abunPlot <- as_tibble(rbind(c(33, "Abundance", "Wildfire", babun.wf$beta, babun.wf$ci.lb, babun.wf$ci.ub, babun.wf$pval, babun.wf$k),
                            c(32, "Abundance", "Prescribed Fire", babun.pf$beta, babun.pf$ci.lb, babun.pf$ci.ub, babun.pf$pval, babun.pf$k),
                            #c(31, "Abundance", "Treatment+Fire", babun.tf$beta, babun.tf$ci.lb, babun.tf$ci.ub, babun.tf$pval, babun.tf$k),
                            c(27, "Abundance", "Growing", babun.g$beta, babun.g$ci.lb, babun.g$ci.ub, babun.g$pval, babun.g$k),
                            c(26, "Abundance", "Dormant", babun.ng$beta, babun.ng$ci.lb, babun.ng$ci.ub, babun.ng$pval, babun.ng$k),
                            c(22, "Abundance", "Broadleaf Forest", babun.df$beta, babun.df$ci.lb, babun.df$ci.ub, babun.df$pval, babun.df$k),
                            c(21, "Abundance", "Coniferous Forest", babun.cf$beta, babun.cf$ci.lb, babun.cf$ci.ub, babun.cf$pval, babun.cf$k),
                            c(20, "Abundance", "Grassland", babun.gr$beta, babun.gr$ci.lb, babun.gr$ci.ub, babun.gr$pval, babun.gr$k),
                            c(16, "Abundance", "Low Severity", babun.lo$beta, babun.lo$ci.lb, babun.lo$ci.ub, babun.lo$pval, babun.lo$k),
                            c(15, "Abundance", "Mid Severity", babun.mi$beta, babun.mi$ci.lb, babun.mi$ci.ub, babun.mi$pval, babun.mi$k),
                            c(14, "Abundance", "High Severity", babun.hi$beta, babun.hi$ci.lb, babun.hi$ci.ub, babun.hi$pval, babun.hi$k),
                            c(10, "Abundance", "0 Years", babun.ze$beta, babun.ze$ci.lb, babun.ze$ci.ub, babun.ze$pval, babun.ze$k),
                            c(9, "Abundance", "1-5 Years", babun.of$beta, babun.of$ci.lb, babun.of$ci.ub, babun.of$pval, babun.of$k),
                            #c(8, "Abundance", "6-10 Years", babun.ft$beta, babun.ft$ci.lb, babun.ft$ci.ub, babun.ft$pval, babun.ft$k),
                            c(7, "Abundance", "11+ Years", babun.tp$beta, babun.tp$ci.lb, babun.tp$ci.ub, babun.tp$pval, babun.tp$k),
                            c(3, "Abundance", "Multiple Fires", babun.ff$beta, babun.ff$ci.lb, babun.ff$ci.ub, babun.ff$pval, babun.ff$k),
                            c(2, "Abundance", "Single Fire", babun.nff$beta, babun.nff$ci.lb, babun.nff$ci.ub, babun.nff$pval, babun.nff$k),
                            c(-2, "Abundance", "Overall Effect", ma01$beta, ma01$ci.lb, ma01$ci.ub, ma01$pval, ma01$k)))
#c(-4, "Abundance", "Pine Barrens", babun.pb$beta, babun.pb$ci.lb, babun.pb$ci.ub, babun.pb$pval, babun.pb$k),
#c(-5, "Abundance", "Valley Forge", babun.vf$beta, babun.vf$ci.lb, babun.vf$ci.ub, babun.vf$pval, babun.vf$k)))

plots.list[[38]] <- ggplot(data=abunPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
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
    ylim(-4,34)+
    scale_x_continuous(limits=c(-5,5), breaks=c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2))+
    ggtitle(label="Abundance", subtitle=paste(dat[1,]$Taxon))


# Richness
richPlot <- as_tibble(rbind(c(33, "Richness", "Wildfire", brich.wf$beta, brich.wf$ci.lb, brich.wf$ci.ub, brich.wf$pval, brich.wf$k),
                            c(32, "Richness", "Prescribed Fire", brich.pf$beta, brich.pf$ci.lb, brich.pf$ci.ub, brich.pf$pval, brich.pf$k),
                            #c(31, "Richness", "Treatment+Fire", brich.tf$beta, brich.tf$ci.lb, brich.tf$ci.ub, brich.tf$pval, brich.tf$k),
                            c(27, "Richness", "Growing", brich.g$beta, brich.g$ci.lb, brich.g$ci.ub, brich.g$pval, brich.g$k),
                            c(26, "Richness", "Dormant", brich.ng$beta, -2, brich.ng$ci.ub, brich.ng$pval, brich.ng$k),
                            c(22, "Richness", "Broadleaf Forest", brich.df$beta, brich.df$ci.lb, brich.df$ci.ub, brich.df$pval, brich.df$k),
                            c(21, "Richness", "Coniferous Forest", brich.cf$beta, brich.cf$ci.lb, brich.cf$ci.ub, brich.cf$pval, brich.cf$k),
                            c(20, "Richness", "Grassland", brich.gr$beta, brich.gr$ci.lb, brich.gr$ci.ub, brich.gr$pval, brich.gr$k),
                            c(16, "Richness", "Low Severity", brich.lo$beta, brich.lo$ci.lb, brich.lo$ci.ub, brich.lo$pval, brich.lo$k),
                            c(15, "Richness", "Mid Severity", brich.mi$beta, brich.mi$ci.lb, brich.mi$ci.ub, brich.mi$pval, brich.mi$k),
                            c(14, "Richness", "High Severity", brich.hi$beta, brich.hi$ci.lb, brich.hi$ci.ub, brich.hi$pval, brich.hi$k),
                            c(10, "Richness", "0 Years", brich.ze$beta, brich.ze$ci.lb, brich.ze$ci.ub, brich.ze$pval, brich.ze$k),
                            c(9, "Richness", "1-5 Years", brich.of$beta, brich.of$ci.lb, brich.of$ci.ub, brich.of$pval, brich.of$k),
                            c(8, "Richness", "6-10 Years", brich.ft$beta, brich.ft$ci.lb, brich.ft$ci.ub, brich.ft$pval, brich.ft$k),
                            #c(7, "Richness", "11+ Years", brich.tp$beta, brich.tp$ci.lb, brich.tp$ci.ub, brich.tp$pval, brich.tp$k),
                            c(3, "Richness", "Multiple Fires", brich.ff$beta, brich.ff$ci.lb, brich.ff$ci.ub, brich.ff$pval, brich.ff$k),
                            c(2, "Richness", "Single Fire", brich.nff$beta, brich.nff$ci.lb, brich.nff$ci.ub, brich.nff$pval, brich.nff$k),
                            c(-2, "Richness", "Overall Effect", ma02$beta, ma02$ci.lb, ma02$ci.ub, ma02$pval, ma02$k)))
#c(-4, "Richness", "Pine Barrens", brich.pb$beta, brich.pb$ci.lb, brich.pb$ci.ub, brich.pb$pval, brich.pb$k),
#c(-5, "Richness", "Valley Forge", brich.vf$beta, brich.vf$ci.lb, brich.vf$ci.ub, brich.vf$pval, brich.vf$k)))

plots.list[[39]] <- ggplot(data=richPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
    geom_vline(xintercept=0, color="gray25", lty="dotted")+
    # geom_segment(data=richPlot,
    #              aes(x=-1.0, xend=-2.0, y=26, yend=26),
    #               arrow=arrow(length=unit(2, "mm")))+
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
    ylim(-4,34)+
    scale_x_continuous(limits=c(-5,5), breaks=c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2))+
    ggtitle(label="Richness", subtitle=paste(dat[1,]$Taxon))

# Print all to PDF
final.plots <- lapply(plyr::compact(plots.list[-c(2)]), FUN=ggplotGrob)

printPath = paste(dat[1,]$Taxon, ".pdf", sep="")
ggsave(printPath, marrangeGrob(grobs=final.plots, nrow=1, ncol=1), 
       width=10, height=8, units="in")

#### Q-STATISTICS ####
babun.growQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", mods=~Growing.non, random=~1|PublicationCode/EffectCode, slab=EffectCode)
brich.growQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", mods=~Growing.non, random=~1|PublicationCode/EffectCode, slab=EffectCode)

babun.typeQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", mods=~FireType, random=~1|PublicationCode/EffectCode, slab=EffectCode)
brich.typeQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", mods=~FireType, random=~1|PublicationCode/EffectCode, slab=EffectCode)

babun.SeverityQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=filter(babun, Severity != "Mid"), method="REML", mods=~Severity, random=~1|PublicationCode/EffectCode, slab=EffectCode)
brich.SeverityQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=filter(brich, Severity != "Mid"), method="REML", mods=~Severity, random=~1|PublicationCode/EffectCode, slab=EffectCode)

babun.yearsQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=filter(babun, TimeClass %in% c("0", "1to5")), method="REML", mods=~TimeClass, random=~1|PublicationCode/EffectCode, slab=EffectCode)
brich.yearsQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=filter(brich, TimeClass %in% c("0", "1to5")), method="REML", mods=~TimeClass, random=~1|PublicationCode/EffectCode, slab=EffectCode)

babun.freqQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", mods=~FireFreq, random=~1|PublicationCode/EffectCode, slab=EffectCode)
brich.freqQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", mods=~FireFreq, random=~1|PublicationCode/EffectCode, slab=EffectCode)

babun.habitatQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=babun, method="REML", mods=~HabitatType, random=~1|PublicationCode/EffectCode, slab=EffectCode)
brich.habitatQ <- rma.mv(yi=EffectSize, V=EffectSizeStandardError, data=brich, method="REML", mods=~HabitatType, random=~1|PublicationCode/EffectCode, slab=EffectCode)

babun.growQ$QM
babun.growQ$QMp
brich.growQ$QM
brich.growQ$QMp

babun.typeQ$QM
babun.typeQ$QMp
brich.typeQ$QM
brich.typeQ$QMp

babun.SeverityQ$QM
babun.SeverityQ$QMp
brich.SeverityQ$QM
brich.SeverityQ$QMp

babun.yearsQ$QM
babun.yearsQ$QMp
brich.yearsQ$QM
brich.yearsQ$QMp

babun.freqQ$QM
babun.freqQ$QMp
brich.freqQ$QM
brich.freqQ$QMp

babun.habitatQ$QM
babun.habitatQ$QMp
brich.habitatQ$QM
brich.habitatQ$QMp

