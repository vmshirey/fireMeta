## This script reads in effect size and standard error information from
## studies on insect responses to fire/other management practices,
## performs 
##
## Author: Vaughn Shirey
##

# load required packages, establish re-utilized variables, custom functions
require(tidyverse); require(bayesmeta); library(meta)

jackknife <- function(effects){
  eff.copy <- data.frame()
  
  dropVec <- unique(effects[,"EffectCode"]) # publication labels
  sens <- list() # list for critical values
  for(i in 1:length(dropVec)){
    eff.copy <- effects %>% filter(EffectCode != dropVec[i])# drop random study
    
    # run meta on reduced dataset
    ma <- rma.mv(yi=EffectSize, 
                 V=eff.copy$EffectSizeStandardError,
                 random=list(~1|PublicationCode,
                             ~1|EffectCode),
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
dat <- read.csv("Bees/beesFinal.csv")

# subset data by taxon/diversity metric
babun <- dat %>% filter(Taxon == "Bees", Metric=="Abundance") # bee abundance
brich <- dat %>% filter(Taxon == "Bees", Metric=="Richness") # bee richness

# run the analysis for abundance
ma01 <- rma.mv(yi=EffectSize, 
               V=babun$EffectSizeStandardError,
               random=list(~1|PublicationCode,
                           ~1|EffectCode),
               tdist=TRUE,
               data=babun, 
               method="REML", 
               slab=EffectCode)

# jackknife analysis for abundance
babun.sens <- jackknife(babun)
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

# plot mean variation
boxplot(estimate)
abline(h=ma01$b, col="red", lwd=2, lty=2)

# plot means, CIs for jackknife and actual effect
qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=ma01$b, color="red")+
  geom_hline(yintercept=ma01$ci.lb, linetype="dashed")+
  geom_hline(yintercept=ma01$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# run the analysis for richness
ma02 <- rma.mv(yi=EffectSize, 
               V=brich$EffectSizeStandardError,
               random=list(~1|PublicationCode,
                           ~1|EffectCode),
               tdist=TRUE,
               data=brich, 
               method="REML", 
               slab=EffectCode)

# jackknife analysis for abundance
brich.sens <- jackknife(brich)
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

# plot mean variation
boxplot(estimate)
abline(h=ma02$b, col="red", lwd=2, lty=2)

# plot means, CIs for jackknife and actual effect
qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=ma02$b, color="red")+
  geom_hline(yintercept=ma02$ci.lb, linetype="dashed")+
  geom_hline(yintercept=ma02$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Fire Severity on Abundance
babun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(Intensity=="High"), slab=EffectCode)
babun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(Intensity=="Mid"), slab=EffectCode)
babun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(Intensity=="Low"), slab=EffectCode)

babun.sens <- jackknife(subset(babun, Intensity=="High"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.hi$b, color="red")+
  geom_hline(yintercept=babun.hi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.hi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")
                   
babun.sens <- jackknife(subset(babun, Intensity=="Mid"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.mi$b, color="red")+
  geom_hline(yintercept=babun.mi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.mi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")
                   
babun.sens <- jackknife(subset(babun, Intensity=="Low"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.lo$b, color="red")+
  geom_hline(yintercept=babun.lo$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.lo$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")


# Fire Severity on Richness
brich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(Intensity=="High"), slab=EffectCode)
brich.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(Intensity=="Mid"), slab=EffectCode)
brich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(Intensity=="Low"), slab=EffectCode)

brich.sens <- jackknife(subset(brich, Intensity=="High"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.hi$b, color="red")+
  geom_hline(yintercept=brich.hi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.hi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, Intensity=="Mid"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.mi$b, color="red")+
  geom_hline(yintercept=brich.mi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.mi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, Intensity=="Low"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.lo$b, color="red")+
  geom_hline(yintercept=brich.lo$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.lo$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Growing vs. Non-growing On Abundance
babun.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(Growing.non=="Growing"), slab=EffectCode)
babun.ng <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(Growing.non=="Dormant"), slab=EffectCode)

babun.sens <- jackknife(subset(babun, Growing.non=="Growing"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.g$b, color="red")+
  geom_hline(yintercept=babun.g$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.g$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, Growing.non=="Dormant"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.ng$b, color="red")+
  geom_hline(yintercept=babun.ng$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.ng$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Growing vs. Non-growing On Richness
brich.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(Growing.non=="Growing"), slab=EffectCode)
brich.ng <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(Growing.non=="Dormant"), slab=EffectCode)

brich.sens <- jackknife(subset(brich, Growning.non=="Growing"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.g$b, color="red")+
  geom_hline(yintercept=brich.g$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.g$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, Growing.non=="Dormant"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.ng$b, color="red")+
  geom_hline(yintercept=brich.ng$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.ng$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Prescribed vs. Wildfire vs. Treatment+Fire on Abundance
babun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(FireType=="Prescribed Fire"), slab=EffectCode)
babun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(FireType=="Wildfire"), slab=EffectCode)
babun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(FireType=="Treatment+Fire"), slab=EffectCode)

babun.sens <- jackknife(subset(babun, FireType=="Prescribed Fire"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.hi$b, color="red")+
  geom_hline(yintercept=babun.hi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.hi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, FireType=="Wildfire"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.mi$b, color="red")+
  geom_hline(yintercept=babun.mi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.mi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, FireType=="Treatment+Fire"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.lo$b, color="red")+
  geom_hline(yintercept=babun.lo$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.lo$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Prescribed vs. Wildfire vs. Treatment+Fire on Richness
brich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(FireType=="Prescribed Fire"), slab=EffectCode)
brich.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(FireType=="Wildfire"), slab=EffectCode)
brich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(FireType=="Treatment+Fire"), slab=EffectCode)

brich.sens <- jackknife(subset(brich, FireType=="Prescribed Fire"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.hi$b, color="red")+
  geom_hline(yintercept=brich.hi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.hi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, FireType=="Wildfire"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.mi$b, color="red")+
  geom_hline(yintercept=brich.mi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.mi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, FireType=="Treatment+Fire"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.lo$b, color="red")+
  geom_hline(yintercept=brich.lo$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.lo$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Time since fire (0, 1to5, 5to10, 10+) on Abundance
babun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(TimeClass=="0"), slab=EffectCode)
babun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(TimeClass=="1to5"), slab=EffectCode)
babun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(TimeClass=="5to10"), slab=EffectCode)
babun.un <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(TimeClass=="10+"), slab=EffectCode)

babun.sens <- jackknife(subset(babun, TimeClass=="0"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.hi$b, color="red")+
  geom_hline(yintercept=babun.hi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.hi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, TimeClass=="1to5"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.mi$b, color="red")+
  geom_hline(yintercept=babun.mi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.mi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, TimeClass=="5to10"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.lo$b, color="red")+
  geom_hline(yintercept=babun.lo$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.lo$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, TimeClass=="10+"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.un$b, color="red")+
  geom_hline(yintercept=babun.un$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.un$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Time since fire (0, 1to5, 5to10, 10+) on Richness
brich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(TimeClass=="0"), slab=EffectCode)
brich.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(TimeClass=="1to5"), slab=EffectCode)
brich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(TimeClass=="5to10"), slab=EffectCode)
brich.un <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(TimeClass=="10+"), slab=EffectCode)

brich.sens <- jackknife(subset(brich, TimeClass=="0"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.hi$b, color="red")+
  geom_hline(yintercept=brich.hi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.hi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, TimeClass=="1to5"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.mi$b, color="red")+
  geom_hline(yintercept=brich.mi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.mi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, TimeClass=="5to10"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.lo$b, color="red")+
  geom_hline(yintercept=brich.lo$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.lo$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, TimeClass=="10+"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.un$b, color="red")+
  geom_hline(yintercept=brich.un$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.un$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Fire Freq (yes/no) On Abundance
babun.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(FireFreq=="Yes"), slab=EffectCode)
babun.ng <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(FireFreq=="No"), slab=EffectCode)

babun.sens <- jackknife(subset(babun, FireFreq=="Yes"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.g$b, color="red")+
  geom_hline(yintercept=babun.g$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.g$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, FireFreq=="No"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.ng$b, color="red")+
  geom_hline(yintercept=babun.ng$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.ng$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Fire Freq (yes/no) On Richness
brich.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(FireFreq=="Yes"), slab=EffectCode)
brich.ng <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(FireFreq=="No"), slab=EffectCode)

brich.sens <- jackknife(subset(brich, Growning.non=="Yes"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.g$b, color="red")+
  geom_hline(yintercept=brich.g$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.g$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, FireFreq=="No"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.ng$b, color="red")+
  geom_hline(yintercept=brich.ng$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.ng$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Habitat (Deciduous Forest, Grassland, Coniferous Forest) on Abundance
babun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(HabitatType=="Deciduous Forest"), slab=EffectCode)
babun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=EffectCode)
babun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=babun, method="SJ", subset=(HabitatType=="Grassland"), slab=EffectCode)

babun.sens <- jackknife(subset(babun, HabitatType=="Deciduous Forest"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.hi$b, color="red")+
  geom_hline(yintercept=babun.hi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.hi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, HabitatType=="Coniferous Forest"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.mi$b, color="red")+
  geom_hline(yintercept=babun.mi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.mi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

babun.sens <- jackknife(subset(babun, HabitatType=="Grassland"))
estimate <- sapply(babun.sens, "[[", 2)
estimate.lb <- sapply(babun.sens, "[[", 4)
estimate.ub <- sapply(babun.sens, "[[", 5)
pubs.out <- sapply(babun.sens, "[[", 7)

qplot(x=seq(1,length(babun.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(babun.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=babun.lo$b, color="red")+
  geom_hline(yintercept=babun.lo$ci.lb, linetype="dashed")+
  geom_hline(yintercept=babun.lo$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

# Habitat (Deciduous Forest, Grassland, Coniferous Forest) on Richness
brich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(HabitatType=="Deciduous Forest"), slab=EffectCode)
brich.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=EffectCode)
brich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=brich, method="SJ", subset=(HabitatType=="Grassland"), slab=EffectCode)

brich.sens <- jackknife(subset(brich, HabitatType=="Deciduous Forest"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.hi$b, color="red")+
  geom_hline(yintercept=brich.hi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.hi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, HabitatType=="Coniferous Forest"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.mi$b, color="red")+
  geom_hline(yintercept=brich.mi$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.mi$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")

brich.sens <- jackknife(subset(brich, HabitatType=="Grassland"))
estimate <- sapply(brich.sens, "[[", 2)
estimate.lb <- sapply(brich.sens, "[[", 4)
estimate.ub <- sapply(brich.sens, "[[", 5)
pubs.out <- sapply(brich.sens, "[[", 7)

qplot(x=seq(1,length(brich.sens)), y=estimate)+
  geom_errorbar(aes(x=seq(1,length(brich.sens)), ymin=estimate.lb, ymax=estimate.ub, width=0.25))+
  geom_hline(yintercept=brich.lo$b, color="red")+
  geom_hline(yintercept=brich.lo$ci.lb, linetype="dashed")+
  geom_hline(yintercept=brich.lo$ci.ub, linetype="dashed")+
  geom_hline(yintercept=mean(estimate), color="blue")