####################################################
## Analysis code for pyroentomology meta-analysis ##
## Vaughn Shirey, 2019                            ##
####################################################

# load libraries
library(meta)
library(metafor)
library(metaviz)
library(forestplot)
library(glmm)
library(tidyverse)

# Read in data from .csv files
carab <- read.csv("Carabidae/carabidData.csv", header=TRUE, sep=",")
bees <- read.csv("Bees/beeData.csv", header=TRUE, sep=",")
leps <- read.csv("Leps/lepData.csv", header=TRUE, sep=",")

# Modify data to conform to particular research question
# carab <- carab[which(carab$YearSinceFire < 2),]
bees <- bees[which(bees$Author != "Simanonok (Thesis)"),]

carab.abun <- carab[carab$Metric == "Abundance",]
carab.rich <- carab[carab$Metric == "Richness",]

bees.abun <- bees[bees$Metric == "Abundance",]
bees.rich <- bees[bees$Metric == "Richness",]

leps.abun <- leps[leps$Metric == "Abundance",]
leps.rich <- leps[leps$Metric == "Richness",]

# Random effects models
carab.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", slab=paste(carab.abun$Author, carab.abun$YearPublished))
carab.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", slab=paste(carab.rich$Author, carab.rich$YearPublished))

bees.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", slab=paste(bees.abun$Author, bees.abun$YearPublished))
bees.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", slab=paste(bees.rich$Author, bees.rich$YearPublished))

leps.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", slab=paste(leps.abun$Author, leps.abun$YearPublished))
leps.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", slab=paste(leps.rich$Author, leps.rich$YearPublished))

# meta regressions for variable effects
carab.abun.mreg <- rma(yi=EffectSize, 
                      sei=EffectSizeStandardError, 
                      mods=~FireType+Intensity+HabitatType+Growing+TimeClass,
                      data=filter(carab.abun, !Intensity=="N/A" & !Growing=="N/A"),
                      test="knha",
                      method="SJ")
carab.abun.mreg

carab.rich.mreg <- rma(yi=EffectSize, 
                      sei=EffectSizeStandardError, 
                      mods=~FireType+Intensity+HabitatType+Growing+TimeClass,
                      data=filter(carab.rich, !Intensity=="N/A" & !Growing=="N/A"),
                      test="knha",
                      method="SJ")
carab.rich.mreg

bees.abun.mreg <- rma(yi=EffectSize, 
                     sei=EffectSizeStandardError, 
                     mods=~FireType+Intensity+HabitatType+Growing+TimeClass,
                     data=filter(bees.abun, !Intensity=="N/A" & !Growing=="N/A"),
                     test="knha",
                     method="SJ")
bees.abun.mreg

bees.rich.mreg <- rma(yi=EffectSize, 
                      sei=EffectSizeStandardError, 
                      mods=~FireType+Intensity+HabitatType+Growing+TimeClass,
                      data=filter(bees.rich, !Intensity=="N/A" & !Growing=="N/A"),
                      test="knha",
                      method="SJ")
bees.rich.mreg

leps.abun.mreg <- rma(yi=EffectSize, 
                      sei=EffectSizeStandardError, 
                      mods=~FireType+Intensity+HabitatType+Growing+TimeClass,
                      data=filter(leps.abun, !Intensity=="N/A" & !Growing=="N/A"),
                      test="knha",
                      method="SJ")
leps.abun.mreg

leps.rich.mreg <- rma(yi=EffectSize, 
                      sei=EffectSizeStandardError, 
                      mods=~FireType+Intensity+HabitatType+Growing+TimeClass,
                      data=filter(leps.rich, !Intensity=="N/A" & !Growing=="N/A"),
                      test="knha",
                      method="SJ")
leps.rich.mreg

# Forest plots for each metric
forest(carab.abun.meta, slab=paste(carab.abun$Author, carab.abun$YearPublished), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)
forest(carab.rich.meta, slab=paste(carab.rich$Author, carab.rich$YearPublished), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)

forest(bees.abun.meta, slab=paste(bees.abun$Author, bees.abun$Year.Published), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = TRUE)
forest(bees.rich.meta, slab=paste(bees.rich$Author, bees.rich$Year.Published), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = TRUE)

forest(leps.abun.meta, slab=paste(leps.abun$Author, leps.abun$Year.Published), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = TRUE)
forest(leps.rich.meta, slab=paste(leps.rich$Author, leps.rich$Year.Published), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = TRUE)

# Subgroup forest plots for each metric (code from metafor website)

# CARABIDS #
carab.abun.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(FireType=="Wildfire"), slab=Author)
carab.abun.pf <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(FireType=="Prescribed Fire"), slab=Author)

carab.rich.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(FireType=="Wildfire"), slab=Author)
carab.rich.pf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(FireType=="Prescribed Fire"), slab=Author)

carab.abun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(Intensity=="High"), slab=Author)
carab.abun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(Intensity=="Mid"), slab=Author)
carab.abun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(Intensity=="Low"), slab=Author)

carab.rich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(Intensity=="High"), slab=Author)
carab.rich.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(Intensity=="Mid"), slab=Author)
carab.rich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(Intensity=="Low"), slab=Author)

carab.abun.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(Growing=="Growing"), slab=Author)
carab.abun.ng <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(Growing=="Dormant"), slab=Author)

carab.rich.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(Growing=="Growing"), slab=Author)
carab.rich.ng <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(Growing=="Dormant"), slab=Author)

carab.abun.df <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(HabitatType=="Deciduous Forest"), slab=Author)
carab.abun.gr <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(HabitatType=="Grassland"), slab=Author)
carab.abun.cf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=Author)
carab.abun.de <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(HabitatType=="Shrubland"), slab=Author)
carab.abun.mf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(HabitatType=="Mixed Forest"), slab=Author)

carab.rich.df <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(HabitatType=="Deciduous Forest"), slab=Author)
carab.rich.gr <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(HabitatType=="Grassland"), slab=Author)
carab.rich.cf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=Author)
carab.rich.de <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(HabitatType=="Shrubland"), slab=Author)
carab.rich.mf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(HabitatType=="Mixed Forest"), slab=Author)

carab.abun.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(TimeClass=="0"), slab=Author)
carab.abun.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(TimeClass=="1-5"), slab=Author)
carab.abun.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(TimeClass=="5+"), slab=Author)

carab.rich.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(TimeClass=="0"), slab=Author)
carab.rich.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(TimeClass=="1-5"), slab=Author)
carab.rich.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(TimeClass=="5+"), slab=Author)

# BEES #
bees.abun.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(FireType=="Wildfire"), slab=Author)
bees.abun.pf <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(FireType=="Prescribed Fire"), slab=Author)

bees.rich.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(FireType=="Wildfire"), slab=Author)
bees.rich.pf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(FireType=="Prescribed Fire"), slab=Author)

bees.abun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Intensity=="High"), slab=Author)
bees.abun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Intensity=="Mid"), slab=Author)
bees.abun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Intensity=="Low"), slab=Author)

bees.rich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Intensity=="High"), slab=Author)
bees.rich.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Intensity=="Mid"), slab=Author)
bees.rich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Intensity=="Low"), slab=Author)

bees.abun.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Growing=="Growing"), slab=Author)
bees.abun.ng <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Growing=="Dormant"), slab=Author)

bees.rich.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Growing=="Growing"), slab=Author)
bees.rich.ng <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Growing=="Dormant"), slab=Author)

bees.abun.df <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(HabitatType=="Deciduous Forest"), slab=Author)
bees.abun.gr <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(HabitatType=="Grassland"), slab=Author)
bees.abun.cf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=Author)
bees.abun.de <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(HabitatType=="Desert"), slab=Author)

bees.rich.df <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(HabitatType=="Deciduous Forest"), slab=Author)
bees.rich.gr <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(HabitatType=="Grassland"), slab=Author)
bees.rich.cf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=Author)
bees.rich.de <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(HabitatType=="Desert"), slab=Author)

bees.abun.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(TimeClass=="0"), slab=Author)
bees.abun.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(TimeClass=="1-5"), slab=Author)
bees.abun.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(TimeClass=="5+"), slab=Author)

bees.rich.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(TimeClass=="0"), slab=Author)
bees.rich.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(TimeClass=="1-5"), slab=Author)
bees.rich.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(TimeClass=="5+"), slab=Author)

bees.abun.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(YearSinceFire <= 3), slab=Author)
bees.abun.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(YearSinceFire > 3), slab=Author)

bees.rich.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(YearSinceFire <= 3), slab=Author)
bees.rich.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(YearSinceFire > 3), slab=Author)

# Subgroup Forest Plots #

bees.richPlot <- as_tibble(rbind(c(19, "Bee Richness", "Wildfire", bees.rich.wf$beta, bees.rich.wf$ci.lb, bees.rich.wf$ci.ub, bees.rich.wf$pval, bees.rich.wf$k),
                             c(18, "Bee Richness", "Prescribed Fire", bees.rich.pf$beta, bees.rich.pf$ci.lb, bees.rich.pf$ci.ub, bees.rich.pf$pval, bees.rich.pf$k),
                             c(15, "Bee Richness", "Growing", bees.rich.g$beta, bees.rich.g$ci.lb, bees.rich.g$ci.ub, bees.rich.g$pval, bees.rich.g$k),
                             c(14, "Bee Richness", "Non-growing", bees.rich.ng$beta, bees.rich.ng$ci.lb, bees.rich.ng$ci.ub, bees.rich.ng$pval, bees.rich.ng$k),
                             c(11, "Bee Richness", "Deciduous Forest", bees.rich.df$beta, bees.rich.df$ci.lb, bees.rich.df$ci.ub, bees.rich.df$pval, bees.rich.df$k),
                             c(10, "Bee Richness", "Coniferous Forest", bees.rich.cf$beta, bees.rich.cf$ci.lb, bees.rich.cf$ci.ub, bees.rich.cf$pval, bees.rich.cf$k),
                             c(9, "Bee Richness", "Grassland", bees.rich.gr$beta, bees.rich.gr$ci.lb, bees.rich.gr$ci.ub, bees.rich.gr$pval, bees.rich.gr$k),
                             c(8, "Bee Richness", "Desert", bees.rich.de$beta, bees.rich.de$ci.lb, bees.rich.de$ci.ub, bees.rich.de$pval, bees.rich.de$k),
                             c(5, "Bee Richness", "Low Severity", bees.rich.lo$beta, bees.rich.lo$ci.lb, bees.rich.lo$ci.ub, bees.rich.lo$pval, bees.rich.lo$k),
                             c(2, "Bee Richness", "0 Years", bees.rich.in$beta, bees.rich.in$ci.lb, bees.rich.in$ci.ub, bees.rich.in$pval, bees.rich.in$k),
                             c(1, "Bee Richness", "1-5 Years", bees.rich.nd$beta, bees.rich.nd$ci.lb, bees.rich.nd$ci.ub, bees.rich.nd$pval, bees.rich.nd$k),
                             c(0, "Bee Richness", "5+ Years", bees.rich.ol$beta, bees.rich.ol$ci.lb, bees.rich.ol$ci.ub, bees.rich.ol$pval, bees.rich.ol$k),
                             c(-2, "Bee Richness", "Overall Effect", bees.rich.meta$beta, bees.rich.meta$ci.lb, bees.rich.meta$ci.ub, bees.rich.meta$pval, bees.rich.meta$k)))

ggplot(data=bees.richPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
  geom_vline(xintercept=0, color="gray25", lty="dotted")+
  geom_errorbarh(aes(xmin=as.numeric(V5), xmax=as.numeric(V6)), width=0.1)+
  geom_point(size=3, pch=22, fill="black")+
  geom_text(aes(x=0, label=paste(V3, " (",V8,")", sep="")), position=position_nudge(x=-3.25))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.line.x = element_line(color="black"))+
  xlab("Effect Size (d)")+
  ylim(-3,21)+
  scale_x_continuous(limits=c(-4,3), breaks=c(-2, -1, -0.5, 0, 0.5, 1, 2))

# Funnel plots for each analysis
funnel(carab.abun.meta)
funnel(carab.rich.meta)

funnel(bees.abun.meta)
funnel(bees.rich.meta)

funnel(leps.abun.meta)
funnel(leps.rich.meta)