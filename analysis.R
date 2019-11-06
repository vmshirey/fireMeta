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
all_taxa <- read.csv("all_taxa/all_taxa.csv", header=TRUE, sep=",")

# Modify data to conform to particular research question
# carab <- carab[which(carab$YearSinceFire < 2),]
bees <- bees[which(bees$Author != "Simanonok (Thesis)"),]

carab.abun <- carab[carab$Metric == "Abundance",]
carab.rich <- carab[carab$Metric == "Richness",]

bees.abun <- bees[bees$Metric == "Abundance",]
bees.rich <- bees[bees$Metric == "Richness",]

leps.abun <- leps[leps$Metric == "Abundance",]
leps.rich <- leps[leps$Metric == "Richness",]

all_taxa.abun <- all_taxa[all_taxa$Metric == "Abundance",]
all_taxa.rich <- all_taxa[all_taxa$Metric == "Richness",]

# Random effects models
carab.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", slab=paste(carab.abun$Author, carab.abun$YearPublished))
carab.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", slab=paste(carab.rich$Author, carab.rich$YearPublished))

bees.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.abun, method="SJ", slab=paste(bees.abun$Author, bees.abun$YearPublished))
bees.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=bees.rich, method="SJ", slab=paste(bees.rich$Author, bees.rich$YearPublished))

leps.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", slab=paste(leps.abun$Author, leps.abun$YearPublished))
leps.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", slab=paste(leps.rich$Author, leps.rich$YearPublished))

all_taxa.abun.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", slab=paste(all_taxa.abun$Author, all_taxa.abun$YearPublished))
all_taxa.rich.meta <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", slab=paste(all_taxa.rich$Author, all_taxa.rich$YearPublished))

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
carab.abun.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(TimeClass=="1to5"), slab=Author)
carab.abun.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.abun, method="SJ", subset=(TimeClass=="5+"), slab=Author)

carab.rich.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(TimeClass=="0"), slab=Author)
carab.rich.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=carab.rich, method="SJ", subset=(TimeClass=="1to5"), slab=Author)
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

# BUTTERFLIES #
leps.abun.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(FireType=="Wildfire"), slab=Author)
leps.abun.pf <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(FireType=="Prescribed"), slab=Author)

leps.rich.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(FireType=="Wildfire"), slab=Author)
leps.rich.pf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(FireType=="Prescribed"), slab=Author)

leps.abun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(Intensity=="High"), slab=Author)
leps.abun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(Intensity=="Mid"), slab=Author)
leps.abun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(Intensity=="Low"), slab=Author)

leps.rich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(Intensity=="High"), slab=Author)
leps.rich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(Intensity=="Low"), slab=Author)

leps.abun.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(Growing=="Growing"), slab=Author)
leps.abun.ng <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(Growing=="Dormant"), slab=Author)

leps.rich.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(Growing=="Growing"), slab=Author)
leps.rich.ng <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(Growing=="Dormant"), slab=Author)

leps.abun.gr <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(HabitatType=="Grassland"), slab=Author)
leps.abun.cf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=Author)

leps.rich.gr <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(HabitatType=="Grassland"), slab=Author)
leps.rich.cf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=Author)

leps.abun.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(TimeClass=="0"), slab=Author)
leps.abun.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(TimeClass=="1to5"), slab=Author)
leps.abun.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.abun, method="SJ", subset=(TimeClass=="5+"), slab=Author)

leps.rich.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(TimeClass=="0"), slab=Author)
leps.rich.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(TimeClass=="1to5"), slab=Author)
leps.rich.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=leps.rich, method="SJ", subset=(TimeClass=="5+"), slab=Author)

## ALL TAXA ##
all_taxa.abun.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(FireType=="Wildfire"), slab=Author)
all_taxa.abun.pf <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(FireType=="Prescribed Fire"), slab=Author)

all_taxa.rich.wf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(FireType=="Wildfire"), slab=Author)
all_taxa.rich.pf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(FireType=="Prescribed Fire"), slab=Author)

all_taxa.abun.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(Intensity=="High"), slab=Author)
all_taxa.abun.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(Intensity=="Mid"), slab=Author)
all_taxa.abun.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(Intensity=="Low"), slab=Author)

all_taxa.rich.hi <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(Intensity=="High"), slab=Author)
all_taxa.rich.mi <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(Intensity=="Mid"), slab=Author)
all_taxa.rich.lo <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(Intensity=="Low"), slab=Author)

all_taxa.abun.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(Growing=="Growing"), slab=Author)
all_taxa.abun.ng <-rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(Growing=="Dormant"), slab=Author)

all_taxa.rich.g <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(Growing=="Growing"), slab=Author)
all_taxa.rich.ng <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(Growing=="Dormant"), slab=Author)

all_taxa.abun.df <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(HabitatType=="Deciduous Forest"), slab=Author)
all_taxa.abun.gr <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(HabitatType=="Grassland"), slab=Author)
all_taxa.abun.cf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=Author)
all_taxa.abun.de <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(HabitatType=="Shrubland"), slab=Author)
all_taxa.abun.mf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(HabitatType=="Mixed Forest"), slab=Author)

all_taxa.rich.df <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(HabitatType=="Deciduous Forest"), slab=Author)
all_taxa.rich.gr <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(HabitatType=="Grassland"), slab=Author)
all_taxa.rich.cf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(HabitatType=="Coniferous Forest"), slab=Author)
all_taxa.rich.de <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(HabitatType=="Shrubland"), slab=Author)
all_taxa.rich.mf <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(HabitatType=="Mixed Forest"), slab=Author)

all_taxa.abun.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(TimeClass=="0"), slab=Author)
all_taxa.abun.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(TimeClass=="1to5"), slab=Author)
all_taxa.abun.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.abun, method="SJ", subset=(TimeClass=="5+"), slab=Author)

all_taxa.rich.in <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(TimeClass=="0"), slab=Author)
all_taxa.rich.nd <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(TimeClass=="1to5"), slab=Author)
all_taxa.rich.ol <- rma(yi=EffectSize, sei=EffectSizeStandardError, data=all_taxa.rich, method="SJ", subset=(TimeClass=="5+"), slab=Author)


# Subgroup Forest Plots #
## BEES RICHNESS PLOT ##
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

## BEES ABUNDANCE PLOT ##
bees.abunPlot <- as_tibble(rbind(c(19, "Bee Abundance", "Wildfire", bees.abun.wf$beta, bees.abun.wf$ci.lb, bees.abun.wf$ci.ub, bees.abun.wf$pval, bees.abun.wf$k),
                                 c(18, "Bee Abundance", "Prescribed Fire", bees.abun.pf$beta, bees.abun.pf$ci.lb, bees.abun.pf$ci.ub, bees.abun.pf$pval, bees.abun.pf$k),
                                 c(15, "Bee Abundance", "Growing", bees.abun.g$beta, bees.abun.g$ci.lb, bees.abun.g$ci.ub, bees.abun.g$pval, bees.abun.g$k),
                                 c(14, "Bee Abundance", "Non-growing", bees.abun.ng$beta, bees.abun.ng$ci.lb, bees.abun.ng$ci.ub, bees.abun.ng$pval, bees.abun.ng$k),
                                 c(11, "Bee Abundance", "Deciduous Forest", bees.abun.df$beta, bees.abun.df$ci.lb, bees.abun.df$ci.ub, bees.abun.df$pval, bees.abun.df$k),
                                 c(10, "Bee Abundance", "Coniferous Forest", bees.abun.cf$beta, bees.abun.cf$ci.lb, bees.abun.cf$ci.ub, bees.abun.cf$pval, bees.abun.cf$k),
                                 c(9, "Bee Abundance", "Grassland", bees.abun.gr$beta, bees.abun.gr$ci.lb, bees.abun.gr$ci.ub, bees.abun.gr$pval, bees.abun.gr$k),
                                 c(5, "Bee Abundance", "Low Severity", bees.abun.lo$beta, bees.abun.lo$ci.lb, bees.abun.lo$ci.ub, bees.abun.lo$pval, bees.abun.lo$k),
                                 c(2, "Bee Abundance", "0 Years", bees.abun.in$beta, bees.abun.in$ci.lb, bees.abun.in$ci.ub, bees.abun.in$pval, bees.abun.in$k),
                                 c(1, "Bee Abundance", "1-5 Years", bees.abun.nd$beta, bees.abun.nd$ci.lb, bees.abun.nd$ci.ub, bees.abun.nd$pval, bees.abun.nd$k),
                                 c(0, "Bee Abundance", "5+ Years", bees.abun.ol$beta, bees.abun.ol$ci.lb, bees.abun.ol$ci.ub, bees.abun.ol$pval, bees.abun.ol$k),
                                 c(-2, "Bee Abundance", "Overall Effect", bees.abun.meta$beta, bees.abun.meta$ci.lb, bees.abun.meta$ci.ub, bees.abun.meta$pval, bees.abun.meta$k)))

ggplot(data=bees.abunPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
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

## CARABID RICHNESS PLOT ##
carab.richPlot <- as_tibble(rbind(c(19, "Carabid Richness", "Wildfire", carab.rich.wf$beta, carab.rich.wf$ci.lb, carab.rich.wf$ci.ub, carab.rich.wf$pval, carab.rich.wf$k),
                                 c(18, "Carabid Richness", "Prescribed Fire", carab.rich.pf$beta, carab.rich.pf$ci.lb, carab.rich.pf$ci.ub, carab.rich.pf$pval, carab.rich.pf$k),
                                 c(15, "Carabid Richness", "Growing", carab.rich.g$beta, carab.rich.g$ci.lb, carab.rich.g$ci.ub, carab.rich.g$pval, carab.rich.g$k),
                                 c(14, "Carabid Richness", "Non-growing", carab.rich.ng$beta, carab.rich.ng$ci.lb, carab.rich.ng$ci.ub, carab.rich.ng$pval, carab.rich.ng$k),
                                 c(11, "Carabid Richness", "Deciduous Forest", carab.rich.df$beta, carab.rich.df$ci.lb, carab.rich.df$ci.ub, carab.rich.df$pval, carab.rich.df$k),
                                 c(10, "Carabid Richness", "Coniferous Forest", carab.rich.cf$beta, carab.rich.cf$ci.lb, carab.rich.cf$ci.ub, carab.rich.cf$pval, carab.rich.cf$k),
                                 c(8, "Carabid Richness", "Grassland", carab.rich.gr$beta, carab.rich.gr$ci.lb, carab.rich.gr$ci.ub, carab.rich.gr$pval, carab.rich.gr$k),
                                 c(7, "Carabid Richness", "Shrubland", carab.rich.de$beta, carab.rich.de$ci.lb, carab.rich.de$ci.ub, carab.rich.de$pval, carab.rich.de$k),
                                 c(4, "Carabid Richness", "Low Severity", carab.rich.lo$beta, carab.rich.lo$ci.lb, carab.rich.lo$ci.ub, carab.rich.lo$pval, carab.rich.lo$k),
                                 c(3, "Carabid Richness", "Mid Severity", carab.rich.mi$beta, carab.rich.mi$ci.lb, carab.rich.mi$ci.ub, carab.rich.mi$pval, carab.rich.mi$k),
                                 c(2, "Carabid Richness", "High Severity", carab.rich.hi$beta, carab.rich.hi$ci.lb, carab.rich.hi$ci.ub, carab.rich.hi$pval, carab.rich.hi$k),
                                 c(-1, "Carabid Richness", "0 Years", carab.rich.in$beta, carab.rich.in$ci.lb, carab.rich.in$ci.ub, carab.rich.in$pval, carab.rich.in$k),
                                 c(-2, "Carabid Richness", "1-5 Years", carab.rich.nd$beta, carab.rich.nd$ci.lb, carab.rich.nd$ci.ub, carab.rich.nd$pval, carab.rich.nd$k),
                                 c(-3, "Carabid Richness", "5+ Years", carab.rich.ol$beta, carab.rich.ol$ci.lb, carab.rich.ol$ci.ub, carab.rich.ol$pval, carab.rich.ol$k),
                                 c(-5, "Carabid Richness", "Overall Effect", carab.rich.meta$beta, carab.rich.meta$ci.lb, carab.rich.meta$ci.ub, carab.rich.meta$pval, carab.rich.meta$k)))

ggplot(data=carab.richPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
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
  ylim(-6,21)+
  scale_x_continuous(limits=c(-4,3), breaks=c(-2, -1, -0.5, 0, 0.5, 1, 2, 3))

## CARABID ABUNDANCE PLOT ##
carab.abunPlot <- as_tibble(rbind(c(19, "Carabid Richness", "Wildfire", carab.abun.wf$beta, carab.abun.wf$ci.lb, carab.abun.wf$ci.ub, carab.abun.wf$pval, carab.abun.wf$k),
                                  c(18, "Carabid Richness", "Prescribed Fire", carab.abun.pf$beta, carab.abun.pf$ci.lb, carab.abun.pf$ci.ub, carab.abun.pf$pval, carab.abun.pf$k),
                                  c(15, "Carabid Richness", "Growing", carab.abun.g$beta, carab.abun.g$ci.lb, carab.abun.g$ci.ub, carab.abun.g$pval, carab.abun.g$k),
                                  c(14, "Carabid Richness", "Non-growing", carab.abun.ng$beta, carab.abun.ng$ci.lb, carab.abun.ng$ci.ub, carab.abun.ng$pval, carab.abun.ng$k),
                                  c(11, "Carabid Richness", "Deciduous Forest", carab.abun.df$beta, carab.abun.df$ci.lb, carab.abun.df$ci.ub, carab.abun.df$pval, carab.abun.df$k),
                                  c(10, "Carabid Richness", "Coniferous Forest", carab.abun.cf$beta, carab.abun.cf$ci.lb, carab.abun.cf$ci.ub, carab.abun.cf$pval, carab.abun.cf$k),
                                  c(8, "Carabid Richness", "Grassland", carab.abun.gr$beta, carab.abun.gr$ci.lb, carab.abun.gr$ci.ub, carab.abun.gr$pval, carab.abun.gr$k),
                                  c(7, "Carabid Richness", "Shrubland", carab.abun.de$beta, carab.abun.de$ci.lb, carab.abun.de$ci.ub, carab.abun.de$pval, carab.abun.de$k),
                                  c(4, "Carabid Richness", "Low Severity", carab.abun.lo$beta, carab.abun.lo$ci.lb, carab.abun.lo$ci.ub, carab.abun.lo$pval, carab.abun.lo$k),
                                  c(3, "Carabid Richness", "Mid Severity", carab.abun.mi$beta, carab.abun.mi$ci.lb, carab.abun.mi$ci.ub, carab.abun.mi$pval, carab.abun.mi$k),
                                  c(2, "Carabid Richness", "High Severity", carab.abun.hi$beta, carab.abun.hi$ci.lb, carab.abun.hi$ci.ub, carab.abun.hi$pval, carab.abun.hi$k),
                                  c(-1, "Carabid Richness", "0 Years", carab.abun.in$beta, carab.abun.in$ci.lb, carab.abun.in$ci.ub, carab.abun.in$pval, carab.abun.in$k),
                                  c(-2, "Carabid Richness", "1-5 Years", carab.abun.nd$beta, carab.abun.nd$ci.lb, carab.abun.nd$ci.ub, carab.abun.nd$pval, carab.abun.nd$k),
                                  c(-3, "Carabid Richness", "5+ Years", carab.abun.ol$beta, carab.abun.ol$ci.lb, carab.abun.ol$ci.ub, carab.abun.ol$pval, carab.abun.ol$k),
                                  c(-5, "Carabid Richness", "Overall Effect", carab.abun.meta$beta, carab.abun.meta$ci.lb, carab.abun.meta$ci.ub, carab.abun.meta$pval, carab.abun.meta$k)))

ggplot(data=carab.abunPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
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
  ylim(-6,21)+
  scale_x_continuous(limits=c(-4,3), breaks=c(-2, -1, -0.5, 0, 0.5, 1, 2, 3))

## BUTTERFLY RICHNESS PLOT ##
leps.richPlot <- as_tibble(rbind(c(19, "Butterfly Richness", "Wildfire", leps.rich.wf$beta, leps.rich.wf$ci.lb, leps.rich.wf$ci.ub, leps.rich.wf$pval, leps.rich.wf$k),
                                 c(18, "Butterfly Richness", "Prescribed Fire", leps.rich.pf$beta, leps.rich.pf$ci.lb, leps.rich.pf$ci.ub, leps.rich.pf$pval, leps.rich.pf$k),
                                 c(15, "Butterfly Richness", "Growing", leps.rich.g$beta, leps.rich.g$ci.lb, leps.rich.g$ci.ub, leps.rich.g$pval, leps.rich.g$k),
                                 c(14, "Butterfly Richness", "Non-growing", leps.rich.ng$beta, leps.rich.ng$ci.lb, leps.rich.ng$ci.ub, leps.rich.ng$pval, leps.rich.ng$k),
                                 c(11, "Butterfly Richness", "Coniferous Forest", leps.rich.cf$beta, leps.rich.cf$ci.lb, leps.rich.cf$ci.ub, leps.rich.cf$pval, leps.rich.cf$k),
                                 c(10, "Butterfly Richness", "Grassland", leps.rich.gr$beta, leps.rich.gr$ci.lb, leps.rich.gr$ci.ub, leps.rich.gr$pval, leps.rich.gr$k),
                                 c(7, "Butterfly Richness", "Low Severity", leps.rich.lo$beta, leps.rich.lo$ci.lb, leps.rich.lo$ci.ub, leps.rich.lo$pval, leps.rich.lo$k),
                                 c(6, "Butterfly Richness", "High Severity", leps.rich.hi$beta, leps.rich.hi$ci.lb, leps.rich.hi$ci.ub, leps.rich.hi$pval, leps.rich.hi$k),
                                 c(3, "Butterfly Richness", "0 Years", leps.rich.in$beta, leps.rich.in$ci.lb, leps.rich.in$ci.ub, leps.rich.in$pval, leps.rich.in$k),
                                 c(2, "Butterfly Richness", "1-5 Years", leps.rich.nd$beta, leps.rich.nd$ci.lb, leps.rich.nd$ci.ub, leps.rich.nd$pval, leps.rich.nd$k),
                                 c(1, "Butterfly Richness", "5+ Years", leps.rich.ol$beta, leps.rich.ol$ci.lb, leps.rich.ol$ci.ub, leps.rich.ol$pval, leps.rich.ol$k),
                                 c(-2, "Butterfly Richness", "Overall Effect", leps.rich.meta$beta, leps.rich.meta$ci.lb, leps.rich.meta$ci.ub, leps.rich.meta$pval, leps.rich.meta$k)))

ggplot(data=leps.richPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
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
  scale_x_continuous(limits=c(-4,3), breaks=c(-2, -1, -0.5, 0, 0.5, 1, 2))+
  annotate("segment", x=-1.5, xend=1.5, y=14, yend=14, arrow=arrow(length=unit(0.02, "npc")))+
  annotate("segment", x=1.5, xend=-1.5, y=14, yend=14, arrow=arrow(length=unit(0.02, "npc")))+
  annotate("segment", x=-1.5, xend=1.5, y=3, yend=3, arrow=arrow(length=unit(0.02, "npc")))+
  annotate("segment", x=1.5, xend=-1.5, y=3, yend=3, arrow=arrow(length=unit(0.02, "npc")))

## BUTTERFLY ABUNDANCE PLOT ##
leps.abunPlot <- as_tibble(rbind(c(19, "Butterfly Richness", "Wildfire", leps.abun.wf$beta, leps.abun.wf$ci.lb, leps.abun.wf$ci.ub, leps.abun.wf$pval, leps.abun.wf$k),
                                 c(18, "Butterfly Richness", "Prescribed Fire", leps.abun.pf$beta, leps.abun.pf$ci.lb, leps.abun.pf$ci.ub, leps.abun.pf$pval, leps.abun.pf$k),
                                 c(15, "Butterfly Richness", "Growing", leps.abun.g$beta, leps.abun.g$ci.lb, leps.abun.g$ci.ub, leps.abun.g$pval, leps.abun.g$k),
                                 c(14, "Butterfly Richness", "Non-growing", leps.abun.ng$beta, leps.abun.ng$ci.lb, leps.abun.ng$ci.ub, leps.abun.ng$pval, leps.abun.ng$k),
                                 c(11, "Butterfly Richness", "Coniferous Forest", leps.abun.cf$beta, leps.abun.cf$ci.lb, leps.abun.cf$ci.ub, leps.abun.cf$pval, leps.abun.cf$k),
                                 c(10, "Butterfly Richness", "Grassland", leps.abun.gr$beta, leps.abun.gr$ci.lb, leps.abun.gr$ci.ub, leps.abun.gr$pval, leps.abun.gr$k),
                                 c(7, "Butterfly Richness", "Low Severity", leps.abun.lo$beta, leps.abun.lo$ci.lb, leps.abun.lo$ci.ub, leps.abun.lo$pval, leps.abun.lo$k),
                                 c(6, "Butterfly Richness", "Mid Severity", leps.abun.mi$beta, leps.abun.mi$ci.lb, leps.abun.mi$ci.ub, leps.abun.mi$pval, leps.abun.mi$k),
                                 c(5, "Butterfly Richness", "High Severity", leps.abun.hi$beta, leps.abun.hi$ci.lb, leps.abun.hi$ci.ub, leps.abun.hi$pval, leps.abun.hi$k),
                                 c(2, "Butterfly Richness", "0 Years", leps.abun.in$beta, leps.abun.in$ci.lb, leps.abun.in$ci.ub, leps.abun.in$pval, leps.abun.in$k),
                                 c(1, "Butterfly Richness", "1-5 Years", leps.abun.nd$beta, leps.abun.nd$ci.lb, leps.abun.nd$ci.ub, leps.abun.nd$pval, leps.abun.nd$k),
                                 c(0, "Butterfly Richness", "5+ Years", leps.abun.ol$beta, leps.abun.ol$ci.lb, leps.abun.ol$ci.ub, leps.abun.ol$pval, leps.abun.ol$k),
                                 c(-3, "Butterfly Richness", "Overall Effect", leps.abun.meta$beta, leps.abun.meta$ci.lb, leps.abun.meta$ci.ub, leps.abun.meta$pval, leps.abun.meta$k)))

ggplot(data=leps.abunPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
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
  scale_x_continuous(limits=c(-4,3), breaks=c(-2, -1, -0.5, 0, 0.5, 1, 2))+
  annotate("segment", x=1.5, xend=-1.5, y=14, yend=14, arrow=arrow(length=unit(0.02, "npc")))+
  annotate("segment", x=-1.5, xend=2.0, y=14, yend=14, arrow=arrow(length=unit(0.02, "npc")))

## COMBINED PLOTS ##
all_taxa.richPlot <- as_tibble(rbind(c(19, "all_taxaid Richness", "Wildfire", all_taxa.rich.wf$beta, all_taxa.rich.wf$ci.lb, all_taxa.rich.wf$ci.ub, all_taxa.rich.wf$pval, all_taxa.rich.wf$k),
                                  c(18, "all_taxaid Richness", "Prescribed Fire", all_taxa.rich.pf$beta, all_taxa.rich.pf$ci.lb, all_taxa.rich.pf$ci.ub, all_taxa.rich.pf$pval, all_taxa.rich.pf$k),
                                  c(15, "all_taxaid Richness", "Growing", all_taxa.rich.g$beta, all_taxa.rich.g$ci.lb, all_taxa.rich.g$ci.ub, all_taxa.rich.g$pval, all_taxa.rich.g$k),
                                  c(14, "all_taxaid Richness", "Non-growing", all_taxa.rich.ng$beta, all_taxa.rich.ng$ci.lb, all_taxa.rich.ng$ci.ub, all_taxa.rich.ng$pval, all_taxa.rich.ng$k),
                                  c(11, "all_taxaid Richness", "Deciduous Forest", all_taxa.rich.df$beta, all_taxa.rich.df$ci.lb, all_taxa.rich.df$ci.ub, all_taxa.rich.df$pval, all_taxa.rich.df$k),
                                  c(10, "all_taxaid Richness", "Coniferous Forest", all_taxa.rich.cf$beta, all_taxa.rich.cf$ci.lb, all_taxa.rich.cf$ci.ub, all_taxa.rich.cf$pval, all_taxa.rich.cf$k),
                                  c(8, "all_taxaid Richness", "Grassland", all_taxa.rich.gr$beta, all_taxa.rich.gr$ci.lb, all_taxa.rich.gr$ci.ub, all_taxa.rich.gr$pval, all_taxa.rich.gr$k),
                                  c(7, "all_taxaid Richness", "Shrubland", all_taxa.rich.de$beta, all_taxa.rich.de$ci.lb, all_taxa.rich.de$ci.ub, all_taxa.rich.de$pval, all_taxa.rich.de$k),
                                  c(4, "all_taxaid Richness", "Low Severity", all_taxa.rich.lo$beta, all_taxa.rich.lo$ci.lb, all_taxa.rich.lo$ci.ub, all_taxa.rich.lo$pval, all_taxa.rich.lo$k),
                                  c(3, "all_taxaid Richness", "Mid Severity", all_taxa.rich.mi$beta, all_taxa.rich.mi$ci.lb, all_taxa.rich.mi$ci.ub, all_taxa.rich.mi$pval, all_taxa.rich.mi$k),
                                  c(2, "all_taxaid Richness", "High Severity", all_taxa.rich.hi$beta, all_taxa.rich.hi$ci.lb, all_taxa.rich.hi$ci.ub, all_taxa.rich.hi$pval, all_taxa.rich.hi$k),
                                  c(-1, "all_taxaid Richness", "0 Years", all_taxa.rich.in$beta, all_taxa.rich.in$ci.lb, all_taxa.rich.in$ci.ub, all_taxa.rich.in$pval, all_taxa.rich.in$k),
                                  c(-2, "all_taxaid Richness", "1-5 Years", all_taxa.rich.nd$beta, all_taxa.rich.nd$ci.lb, all_taxa.rich.nd$ci.ub, all_taxa.rich.nd$pval, all_taxa.rich.nd$k),
                                  c(-3, "all_taxaid Richness", "5+ Years", all_taxa.rich.ol$beta, all_taxa.rich.ol$ci.lb, all_taxa.rich.ol$ci.ub, all_taxa.rich.ol$pval, all_taxa.rich.ol$k),
                                  c(-5, "all_taxaid Richness", "Overall Effect", all_taxa.rich.meta$beta, all_taxa.rich.meta$ci.lb, all_taxa.rich.meta$ci.ub, all_taxa.rich.meta$pval, all_taxa.rich.meta$k)))

ggplot(data=all_taxa.richPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
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
  ylim(-6,21)+
  scale_x_continuous(limits=c(-4,3), breaks=c(-2, -1, -0.5, 0, 0.5, 1, 2, 3))

## all_taxaID ABUNDANCE PLOT ##
all_taxa.abunPlot <- as_tibble(rbind(c(19, "all_taxaid Richness", "Wildfire", all_taxa.abun.wf$beta, all_taxa.abun.wf$ci.lb, all_taxa.abun.wf$ci.ub, all_taxa.abun.wf$pval, all_taxa.abun.wf$k),
                                  c(18, "all_taxaid Richness", "Prescribed Fire", all_taxa.abun.pf$beta, all_taxa.abun.pf$ci.lb, all_taxa.abun.pf$ci.ub, all_taxa.abun.pf$pval, all_taxa.abun.pf$k),
                                  c(15, "all_taxaid Richness", "Growing", all_taxa.abun.g$beta, all_taxa.abun.g$ci.lb, all_taxa.abun.g$ci.ub, all_taxa.abun.g$pval, all_taxa.abun.g$k),
                                  c(14, "all_taxaid Richness", "Non-growing", all_taxa.abun.ng$beta, all_taxa.abun.ng$ci.lb, all_taxa.abun.ng$ci.ub, all_taxa.abun.ng$pval, all_taxa.abun.ng$k),
                                  c(11, "all_taxaid Richness", "Deciduous Forest", all_taxa.abun.df$beta, all_taxa.abun.df$ci.lb, all_taxa.abun.df$ci.ub, all_taxa.abun.df$pval, all_taxa.abun.df$k),
                                  c(10, "all_taxaid Richness", "Coniferous Forest", all_taxa.abun.cf$beta, all_taxa.abun.cf$ci.lb, all_taxa.abun.cf$ci.ub, all_taxa.abun.cf$pval, all_taxa.abun.cf$k),
                                  c(8, "all_taxaid Richness", "Grassland", all_taxa.abun.gr$beta, all_taxa.abun.gr$ci.lb, all_taxa.abun.gr$ci.ub, all_taxa.abun.gr$pval, all_taxa.abun.gr$k),
                                  c(7, "all_taxaid Richness", "Shrubland", all_taxa.abun.de$beta, all_taxa.abun.de$ci.lb, all_taxa.abun.de$ci.ub, all_taxa.abun.de$pval, all_taxa.abun.de$k),
                                  c(4, "all_taxaid Richness", "Low Severity", all_taxa.abun.lo$beta, all_taxa.abun.lo$ci.lb, all_taxa.abun.lo$ci.ub, all_taxa.abun.lo$pval, all_taxa.abun.lo$k),
                                  c(3, "all_taxaid Richness", "Mid Severity", all_taxa.abun.mi$beta, all_taxa.abun.mi$ci.lb, all_taxa.abun.mi$ci.ub, all_taxa.abun.mi$pval, all_taxa.abun.mi$k),
                                  c(2, "all_taxaid Richness", "High Severity", all_taxa.abun.hi$beta, all_taxa.abun.hi$ci.lb, all_taxa.abun.hi$ci.ub, all_taxa.abun.hi$pval, all_taxa.abun.hi$k),
                                  c(-1, "all_taxaid Richness", "0 Years", all_taxa.abun.in$beta, all_taxa.abun.in$ci.lb, all_taxa.abun.in$ci.ub, all_taxa.abun.in$pval, all_taxa.abun.in$k),
                                  c(-2, "all_taxaid Richness", "1-5 Years", all_taxa.abun.nd$beta, all_taxa.abun.nd$ci.lb, all_taxa.abun.nd$ci.ub, all_taxa.abun.nd$pval, all_taxa.abun.nd$k),
                                  c(-3, "all_taxaid Richness", "5+ Years", all_taxa.abun.ol$beta, all_taxa.abun.ol$ci.lb, all_taxa.abun.ol$ci.ub, all_taxa.abun.ol$pval, all_taxa.abun.ol$k),
                                  c(-5, "all_taxaid Richness", "Overall Effect", all_taxa.abun.meta$beta, all_taxa.abun.meta$ci.lb, all_taxa.abun.meta$ci.ub, all_taxa.abun.meta$pval, all_taxa.abun.meta$k)))

ggplot(data=all_taxa.abunPlot, mapping=aes(x=as.numeric(V4), y=as.numeric(V1)))+
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
  ylim(-6,21)+
  scale_x_continuous(limits=c(-4,3), breaks=c(-2, -1, -0.5, 0, 0.5, 1, 2, 3))


## PUBLICATION BIAS ##
# Funnel plots for each analysis
cam.tf <- trimfill(carab.abun.meta)
funnel(cam.tf)

crm.tf <- trimfill(carab.rich.meta)
funnel(crm.tf)

funnel(carab.abun.meta)
funnel(carab.rich.meta)

bam.tf <- trimfill(bees.abun.meta)
funnel(bam.tf)

brm.tf <- trimfill(bees.rich.meta)
funnel(brm.tf)

funnel(bees.abun.meta)
funnel(bees.rich.meta)

lam.tf <- trimfill(leps.abun.meta)
funnel(lam.tf)

lrm.tf <- trimfill(leps.rich.meta)
funnel(lrm.tf)

funnel(leps.abun.meta)
funnel(leps.rich.meta)

atam.tf <- trimfill(all_taxa.abun.meta)
funnel(atam.tf)

atrm.tf <- trimfill(all_taxa.rich.meta)
funnel(atrm.tf)

funnel(all_taxa.abun.meta)
funnel(all_taxa.rich.meta)

# Fail safe number
fsn(bees.abun$EffectSize, sei=bees.abun$EffectSizeStandardError, type="Rosenberg")
fsn(bees.rich$EffectSize, sei=bees.rich$EffectSizeStandardError, type="Rosenberg")

fsn(carab.abun$EffectSize, sei=carab.abun$EffectSizeStandardError, type="Rosenberg")
fsn(carab.rich$EffectSize, sei=carab.rich$EffectSizeStandardError, type="Rosenberg")

fsn(leps.abun$EffectSize, sei=leps.abun$EffectSizeStandardError, type="Rosenberg")
fsn(leps.rich$EffectSize, sei=leps.rich$EffectSizeStandardError, type="Rosenberg")

fsn(all_taxa.abun$EffectSize, sei=all_taxa.abun$EffectSizeStandardError, type="Rosenberg")
fsn(all_taxa.rich$EffectSize, sei=all_taxa.rich$EffectSizeStandardError, type="Rosenberg")
