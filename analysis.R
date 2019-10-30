####################################################
## Analysis code for pyroentomology meta-analysis ##
## Vaughn Shirey, 2019                            ##
####################################################

# load libraries
library(dplyr)
library(meta)
library(metafor)
library(metaviz)
library(forestplot)
library(glmm)

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

# CARABIDS #
par(mfrow=c(1,1), font=1)

plot.data <- structure(list( # wildfire vs prescribed fire
  mean = c(NA, carab.abun.wf$beta, carab.rich.wf$beta, carab.abun.pf$beta, carab.rich.pf$beta),
  lower = c(NA, carab.abun.wf$ci.lb, carab.rich.wf$ci.lb, carab.abun.pf$ci.lb, carab.rich.pf$ci.lb),
  upper = c(NA, carab.abun.wf$ci.ub, carab.rich.wf$ci.ub, carab.abun.pf$ci.ub, carab.rich.pf$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Fire Type", "Wildfire (Abundance) [14]", "Wildfire (Richness) [11]", "Prescribed Fire (Abundance) [15]", "Prescribed Fire (Richness) [13]"),
  c("p-value", format(carab.abun.wf$pval, digits=2), format(carab.rich.wf$pval, digits=2), format(carab.abun.pf$pval, digits=2), format(carab.rich.pf$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # fire intensity
  mean = c(NA, carab.abun.hi$beta, carab.rich.hi$beta, carab.abun.mi$beta, carab.rich.mi$beta, carab.abun.lo$beta, carab.rich.lo$beta),
  lower = c(NA, carab.abun.hi$ci.lb, carab.rich.hi$ci.lb, carab.abun.mi$ci.lb, carab.rich.mi$ci.lb, carab.abun.lo$ci.lb, carab.rich.lo$ci.lb),
  upper = c(NA, carab.abun.hi$ci.ub, carab.rich.hi$ci.ub, carab.abun.mi$ci.ub, carab.rich.mi$ci.ub, carab.abun.lo$ci.ub, carab.rich.lo$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Fire Intensity", "High Intensity (Abundance) [1]", "High Intensity (Richness) [1]", "Mid Intensity (Abundance) [1]", "Mid Intensity (Richness) [1]",
    "Low Intensity (Abundance) [4]", "Low Intensity (Richness) [4]"),
  c("p-value", format(carab.abun.hi$pval, digits=2), format(carab.rich.hi$pval, digits=2), format(carab.abun.mi$pval, digits=2), format(carab.rich.hi$pval, digits=2),
    format(carab.abun.lo$pval, digits=2), format(carab.rich.lo$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # dormant vs. growing season
  mean = c(NA, carab.abun.g$beta, carab.rich.g$beta, carab.abun.ng$beta, carab.rich.ng$beta),
  lower = c(NA, carab.abun.g$ci.lb, carab.rich.g$ci.lb, carab.abun.ng$ci.lb, carab.rich.ng$ci.lb),
  upper = c(NA, carab.abun.g$ci.ub, carab.rich.g$ci.ub, carab.abun.ng$ci.ub, carab.rich.ng$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Season of Fire", "Growing (Abundance) [18]", "Growing (Richness) [15]", "Dormant (Abundance) [1]", "Dormant (Richness) [1]"),
  c("p-value", format(carab.abun.g$pval, digits=2), format(carab.rich.g$pval, digits=2), format(carab.abun.ng$pval, digits=2), format(carab.rich.g$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # habitat type
  mean = c(NA, carab.abun.df$beta, carab.rich.df$beta, carab.abun.cf$beta, carab.rich.cf$beta, carab.abun.mf$beta, carab.rich.mf$beta, carab.abun.gr$beta, carab.rich.gr$beta, carab.abun.de$beta, carab.rich.de$beta),
  lower = c(NA, carab.abun.df$ci.lb, carab.rich.df$ci.lb, carab.abun.cf$ci.lb, carab.rich.cf$ci.lb, carab.abun.mf$ci.lb, carab.rich.mf$ci.lb, carab.abun.gr$ci.lb, carab.rich.gr$ci.lb, carab.abun.de$ci.lb, carab.rich.de$ci.lb),
  upper = c(NA, carab.abun.df$ci.ub, carab.rich.df$ci.ub, carab.abun.cf$ci.ub, carab.rich.cf$ci.ub, carab.abun.mf$ci.ub, carab.rich.mf$ci.ub, carab.abun.gr$ci.ub, carab.rich.gr$ci.ub, carab.abun.de$ci.ub, carab.rich.de$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Habitat of Fire", "Deciduous Forest (Abundance) [2]", "Deciduous Forest (Richness) [2]", "Coniferous Forest (Abundance) [15]", "Coniferous Forest (Richness) [10]",
    "Mixed Forest (Abudance) [1]", "Mixed Forest (Richness) [1]", "Grassland (Abundance) [2]", "Grassland (Richness) [2]",  "Shrubland (Abundance) [2]", "Shrubland (Richness) [2]"),
  c("p-value", format(carab.abun.df$pval, digits=2), format(carab.rich.df$pval, digits=2), format(carab.abun.cf$pval, digits=2), format(carab.rich.cf$pval, digits=2),
    format(carab.abun.mf$pval, digits=2), format(carab.rich.mf$pval, digits=2), format(carab.abun.gr$pval, digits=2), format(carab.rich.gr$pval, digits=2), format(carab.abun.de$pval, digits=2), format(carab.rich.de$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # fire intensity
  mean = c(NA, carab.abun.in$beta, carab.rich.in$beta, carab.abun.nd$beta, carab.rich.nd$beta, carab.abun.ol$beta, carab.rich.ol$beta),
  lower = c(NA, carab.abun.in$ci.lb, carab.rich.in$ci.lb, carab.abun.nd$ci.lb, carab.rich.nd$ci.lb, carab.abun.ol$ci.lb, carab.rich.ol$ci.lb),
  upper = c(NA, carab.abun.in$ci.ub, carab.rich.in$ci.ub, carab.abun.nd$ci.ub, carab.rich.nd$ci.ub, carab.abun.ol$ci.ub, carab.rich.ol$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Time Since Fire", "0 Years (Abundance) [11]", "0 Years (Richness) [9]", "1-5 Years (Abundance) [17]", "1-5 Years (Richness) [14]",
    "5+ Years (Abundance) [1]", "5+ Years (Richness) [1]"),
  c("p-value", format(carab.abun.in$pval, digits=2), format(carab.rich.in$pval, digits=2), format(carab.abun.nd$pval, digits=2), format(carab.rich.nd$pval, digits=2),
    format(carab.abun.ol$pval, digits=2), format(carab.rich.ol$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

# BEES #
par(mfrow=c(1,1), font=1)

plot.data <- structure(list( # wildfire vs prescribed fire
  mean = c(NA, bees.abun.wf$beta, bees.rich.wf$beta, bees.abun.pf$beta, bees.rich.pf$beta),
  lower = c(NA, bees.abun.wf$ci.lb, bees.rich.wf$ci.lb, bees.abun.pf$ci.lb, bees.rich.pf$ci.lb),
  upper = c(NA, bees.abun.wf$ci.ub, bees.rich.wf$ci.ub, bees.abun.pf$ci.ub, bees.rich.pf$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Fire Type", "Wildfire (Abundance) [11]", "Wildfire (Richness) [13]", "Prescribed Fire (Abundance) [5]", "Prescribed Fire (Richness) [3]"),
  c("p-value", format(bees.abun.wf$pval, digits=2), format(bees.rich.wf$pval, digits=2), format(bees.abun.pf$pval, digits=2), format(bees.rich.wf$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # fire intensity
  mean = c(NA, bees.abun.hi$beta, bees.rich.hi$beta, bees.abun.mi$beta, bees.rich.mi$beta, bees.abun.lo$beta, bees.rich.lo$beta),
  lower = c(NA, bees.abun.hi$ci.lb, bees.rich.hi$ci.lb, bees.abun.mi$ci.lb, bees.rich.mi$ci.lb, bees.abun.lo$ci.lb, bees.rich.lo$ci.lb),
  upper = c(NA, bees.abun.hi$ci.ub, bees.rich.hi$ci.ub, bees.abun.mi$ci.ub, bees.rich.mi$ci.ub, bees.abun.lo$ci.ub, bees.rich.lo$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Fire Intensity", "High Intensity (Abundance) [5]", "High Intensity (Richness) [6]", "Mid Intensity (Abundance) [1]", "Mid Intensity (Richness) [4]",
    "Low Intensity (Abundance) [8]", "Low Intensity (Richness) [6]"),
  c("p-value", format(bees.abun.hi$pval, digits=2), format(bees.rich.hi$pval, digits=2), format(bees.abun.mi$pval, digits=2), format(bees.rich.hi$pval, digits=2),
    format(bees.abun.lo$pval, digits=2), format(bees.rich.lo$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # dormant vs. growing season
  mean = c(NA, bees.abun.g$beta, bees.rich.g$beta, bees.abun.ng$beta, bees.rich.ng$beta),
  lower = c(NA, bees.abun.g$ci.lb, bees.rich.g$ci.lb, bees.abun.ng$ci.lb, bees.rich.ng$ci.lb),
  upper = c(NA, bees.abun.g$ci.ub, bees.rich.g$ci.ub, bees.abun.ng$ci.ub, bees.rich.ng$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Season of Fire", "Growing (Abundance) [9]", "Growing (Richness) [8]", "Dormant (Abundance) [5]", "Dormant (Richness) [3]"),
  c("p-value", format(bees.abun.g$pval, digits=2), format(bees.rich.g$pval, digits=2), format(bees.abun.ng$pval, digits=2), format(bees.rich.g$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # habitat type
  mean = c(NA, bees.abun.df$beta, bees.rich.df$beta, bees.abun.cf$beta, bees.rich.cf$beta, bees.abun.gr$beta, bees.rich.gr$beta, bees.rich.de$beta),
  lower = c(NA, bees.abun.df$ci.lb, bees.rich.df$ci.lb, bees.abun.cf$ci.lb, bees.rich.cf$ci.lb, bees.abun.gr$ci.lb, bees.rich.gr$ci.lb, bees.rich.de$ci.lb),
  upper = c(NA, bees.abun.df$ci.ub, bees.rich.df$ci.ub, bees.abun.cf$ci.ub, bees.rich.cf$ci.ub, bees.abun.gr$ci.ub, bees.rich.gr$ci.ub, bees.rich.de$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Habitat of Fire", "Deciduous Forest (Abundance) [7]", "Deciduous Forest (Richness) [5]", "Coniferous Forest (Abundance) [7]", "Coniferous Forest (Richness) [7]",
    "Grassland (Abundance) [2]", "Grassland (Richness) [1]", "Desert (Richness) [3]"),
  c("p-value", format(bees.abun.df$pval, digits=2), format(bees.rich.df$pval, digits=2), format(bees.abun.cf$pval, digits=2), format(bees.rich.cf$pval, digits=2),
    format(bees.abun.gr$pval, digits=2), format(bees.rich.gr$pval, digits=2), format(bees.rich.de$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # fire intensity
  mean = c(NA, bees.abun.in$beta, bees.rich.in$beta, bees.abun.nd$beta, bees.rich.nd$beta, bees.abun.ol$beta, bees.rich.ol$beta),
  lower = c(NA, bees.abun.in$ci.lb, bees.rich.in$ci.lb, bees.abun.nd$ci.lb, bees.rich.nd$ci.lb, bees.abun.ol$ci.lb, bees.rich.ol$ci.lb),
  upper = c(NA, bees.abun.in$ci.ub, bees.rich.in$ci.ub, bees.abun.nd$ci.ub, bees.rich.nd$ci.ub, bees.abun.ol$ci.ub, bees.rich.ol$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Time Since Fire", "0 Years (Abundance) [5]", "0 Years (Richness) [5]", "1-5 Years (Abundance) [4]", "1-5 Years (Richness) [4]",
    "5+ Years (Abundance) [5]", "5+ Years (Richness) [7]"),
  c("p-value", format(bees.abun.in$pval, digits=2), format(bees.rich.in$pval, digits=2), format(bees.abun.nd$pval, digits=2), format(bees.rich.nd$pval, digits=2),
    format(bees.abun.ol$pval, digits=2), format(bees.rich.ol$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

plot.data <- structure(list( # fire intensity
  mean = c(NA, bees.abun.in$beta, bees.rich.in$beta, bees.abun.nd$beta, bees.rich.nd$beta),
  lower = c(NA, bees.abun.in$ci.lb, bees.rich.in$ci.lb, bees.abun.nd$ci.lb, bees.rich.nd$ci.lb),
  upper = c(NA, bees.abun.in$ci.ub, bees.rich.in$ci.ub, bees.abun.nd$ci.ub, bees.rich.nd$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Time Since Fire", "0-5 Years (Abundance) [9]", "0-5 Years (Richness) [9]", "5+ Years (Abundance) [5]", "5+ Years (Richness) [7]"),
  c("p-value", format(bees.abun.in$pval, digits=2), format(bees.rich.in$pval, digits=2), format(bees.abun.nd$pval, digits=2), format(bees.rich.nd$pval, digits=2))
)
forestplot(tabletext, plot.data$mean, plot.data$lower, plot.data$upper, new_page=TRUE)

# Funnel plots for each analysis
funnel(carab.abun.meta)
funnel(carab.rich.meta)

funnel(bees.abun.meta)
funnel(bees.rich.meta)

funnel(leps.abun.meta)
funnel(leps.rich.meta)