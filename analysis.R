####################################################
## Analysis code for pyroentomology meta-analysis ##
## Vaughn Shirey, 2019                            ##
####################################################

# load libraries
library(dplyr)
library(meta)
library(metafor)
library(metaviz)

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

bees.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ")
bees.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ")

# Forest plots for each metric
forest(carab.abun.meta, slab=paste(carab.abun$Author), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)
forest(carab.rich.meta, slab=paste(carab.rich$Author), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = FALSE)

forest(bees.abun.meta, slab=paste(bees.abun$Author), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = TRUE)
forest(bees.rich.meta, slab=paste(bees.rich$Author), order="obs", xlab="Hedge's G", alim=c(-10,10), annotate = TRUE)

# Subgroup forest plots for each metric

bees.abun.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(FireType=="Wildfire"))
bees.abun.pf <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(FireType=="Prescribed Fire"))

bees.rich.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(FireType=="Wildfire"))
bees.rich.pf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(FireType=="Prescribed Fire"))

bees.abun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Intensity=="High"))
bees.abun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Intensity=="Mid"))
bees.abun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Intensity=="Low"))

bees.rich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Intensity=="High"))
bees.rich.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Intensity=="Mid"))
bees.rich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Intensity=="Low"))

bees.abun.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Growing=="Growing"))
bees.abun.ng <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", subset=(Growing=="Dormant"))

bees.rich.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Growing=="Growing"))
bees.rich.ng <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", subset=(Growing=="Dormant"))


# Funnel plots for each analysis
funnel(carab.abun.meta)
funnel(carab.rich.meta)

funnel(bees.abun.meta)
funnel(bees.rich.meta)