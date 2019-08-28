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

# Read in data from .csv files
carab <- read.csv("Carabidae/carabidData.csv", header=TRUE, sep=",")
bees <- read.csv("Bees/beeData.csv", header=TRUE, sep=",")
leps <- NA

# Modify data to conform to particular research question
# carab <- carab[which(carab$YearSinceFire < 2),]
# bees <- bees[which(bees$YearSinceFire < 2),]

carab.abun <- carab[carab$Metric == "Abundance",]
carab.rich <- carab[carab$Metric == "Richness",]

bees.abun <- bees[bees$Metric == "Abundance",]
bees.rich <- bees[bees$Metric == "Richness",]

# Random effects models
carab.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ")
carab.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ")

bees.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", slab=Author)
bees.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", slab=Author)

# Forest plots for each metric
forest(carab.abun.meta, slab=paste(carab.abun$Author), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)
forest(carab.rich.meta, slab=paste(carab.rich$Author), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)

forest(bees.abun.meta, slab=paste(bees.abun$Author), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = TRUE)
forest(bees.rich.meta, slab=paste(bees.rich$Author), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = TRUE)

# Subgroup forest plots for each metric (code from metafor website)

par(mfrow=c(2,1), font=1)

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

# Subgroup Forest Plots #

plot.data <- structure(list( # wildfire vs prescribed fire
  mean = c(NA, bees.abun.wf$beta, bees.rich.wf$beta, bees.abun.pf$beta, bees.rich.pf$beta),
  lower = c(NA, bees.abun.wf$ci.lb, bees.rich.wf$ci.lb, bees.abun.pf$ci.lb, bees.rich.pf$ci.lb),
  upper = c(NA, bees.abun.wf$ci.ub, bees.rich.wf$ci.ub, bees.abun.pf$ci.ub, bees.rich.pf$ci.ub)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -11L),
  class = "data.frame")
tabletext <- cbind(
  c("Fire Type", "Wildfire (Abundance)", "Wildfire (Richness)", "Prescribed Fire (Abundance)", "Prescribed Fire (Richness)"),
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

# Funnel plots for each analysis
funnel(carab.abun.meta)
funnel(carab.rich.meta)

funnel(bees.abun.meta)
funnel(bees.rich.meta)