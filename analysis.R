####################################################
## Analysis code for pyroentomology meta-analysis ##
## Vaughn Shirey, 2019                            ##
####################################################

# load libraries
library(dplyr)
library(meta)
library(metafor)
library(metaviz)

# additional functions
spot.outliers.random<-function(data){
  data<-data
  Author<-data$studlab
  lowerci<-data$lower
  upperci<-data$upper
  m.outliers<-data.frame(Author,lowerci,upperci)
  te.lower<-data$lower.random
  te.upper<-data$upper.random
  dplyr::filter(m.outliers,upperci < te.lower)
  dplyr::filter(m.outliers,lowerci > te.upper)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

influence.analysis<-function(data,method.tau,hakn){
  
  influence.data<-data
  TE<-data$TE
  seTE<-data$seTE
  method.tau<-method.tau
  hakn<-hakn
  
  if(hakn == TRUE){
    res <- rma(yi=TE, sei=seTE, measure="ZCOR", 
               data=influence.data, 
               method = paste(method.tau),
               test="knha")
    res
    inf <- influence(res)
    influence.data<-metainf(data)
    influence.data$I2<-format(round(influence.data$I2,2),nsmall=2)
    plot(inf)
    baujat(data)
    forest(influence.data,
           sortvar=I2,
           rightcols = c("TE","ci","I2"),
           smlab = "Sorted by I-squared")
    forest(influence.data,
           sortvar=TE,
           rightcols = c("TE","ci","I2"),
           smlab = "Sorted by Effect size")
    
  } else {
    
    res <- rma(yi=TE, sei=seTE, measure="ZCOR", 
               data=influence.data, 
               method = paste(method.tau))
    res
    inf <- influence(res)
    influence.data<-metainf(data)
    influence.data$I2<-format(round(influence.data$I2,2),nsmall=2)
    plot(inf)
    baujat(data)
    forest(influence.data,
           sortvar=I2,
           rightcols = c("TE","ci","I2"),
           smlab = "Sorted by I-squared")
    forest(influence.data,
           sortvar=TE,
           rightcols = c("TE","ci","I2"),
           smlab = "Sorted by Effect size")
  }} 

##################################
## Read in data from .csv files ##
##################################

carabidData <- read.csv("Carabidae/carabidData.csv", header=TRUE)
beeData <- read.csv("Bees/beeData.csv", header=TRUE)
butterflyData <- NA

# filter data to only include 0-1 years post fire sampling
carabidData <- carabidData[carabidData$YearSinceFire < 5,]

##########
## BEES ##
##########

# pooled effect size using Hartung-Knapp-Sidik-Jonkman Method #
# abundance
m.abundance.hksj <- metagen(EffectSize, # effect sizes
                            EffectSizeStandardError, # effect size standard errors
                            data = beeData[beeData$Metric=="Abundance",], # data
                            studlab = paste(PaperID), # paper ID for random effect
                            comb.fixed = FALSE, # not doing fixed effects
                            comb.random = TRUE, # doing random effects
                            method.tau = "ML", # Sidik-Jonkman estimator for variance
                            hakn = TRUE, # HKSJ Method
                            prediction = TRUE, # predict an interval
                            sm="SMD") # mean difference effect size

# richness
m.richness.hksj <- metagen(EffectSize, # effect sizes
                           EffectSizeStandardError, # effect size standard errors
                           data = beeData[beeData$Metric=="Richness",], # data
                           studlab = paste(PaperID), # paper ID for random effect
                           comb.fixed = FALSE, # not doing fixed effects
                           comb.random = TRUE, # doing random effects
                           method.tau = "ML", # Sidik-Jonkman estimator for variance
                           hakn = TRUE, # HKSJ Method
                           prediction = TRUE, # predict an interval
                           sm="SMD") # mean difference effect size

# search for outliers and influence #
spot.outliers.random(data=m.abundance.hksj)
spot.outliers.random(data=m.richness.hksj)

## REMOVE THE OUTLIERS ##

# forest plots 

forest(m.abundance.hksj, digits.sd = 2, digits.se = 2)
forest(m.richness.hksj, digits.sd = 2, digits.se = 2)

# subgroup analysis

# abundance
intensity.abundance.subgroup <- update.meta(m.abundance.hksj,
                                            byvar = Intensity,
                                            comb.random = TRUE,
                                            comb.fixed = FALSE)
habitat.abundance.subgroup <- update.meta(m.abundance.hksj,
                                          byvar = HabitatType,
                                          comb.random = TRUE,
                                          comb.fixed = FALSE)
firetype.abundance.subgroup <- update.meta(m.abundance.hksj,
                                           byvar = FireType,
                                           comb.random = TRUE,
                                           comb.fixed = FALSE)
region.abundance.subgroup <- update.meta(m.abundance.hksj,
                                         byvar = Region,
                                         comb.random = TRUE,
                                         comb.fixed = FALSE)
biome.abundance.subgroup <- update.meta(m.abundance.hksj,
                                        byvar = Biome,
                                        comb.random = TRUE,
                                        comb.fixed = FALSE)

forest(intensity.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(habitat.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(firetype.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(region.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(biome.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)

# richness
intensity.richness.subgroup <- update.meta(m.richness.hksj,
                                           byvar = Intensity,
                                           comb.random = TRUE,
                                           comb.fixed = FALSE)
habitat.richness.subgroup <- update.meta(m.richness.hksj,
                                         byvar = HabitatType,
                                         comb.random = TRUE,
                                         comb.fixed = FALSE)
firetype.richness.subgroup <- update.meta(m.richness.hksj,
                                          byvar = FireType,
                                          comb.random = TRUE,
                                          comb.fixed = FALSE)
region.richness.subgroup <- update.meta(m.richness.hksj,
                                        byvar = Region,
                                        comb.random = TRUE,
                                        comb.fixed = FALSE)
biome.richness.subgroup <- update.meta(m.richness.hksj,
                                       byvar = Biome,
                                       comb.random = TRUE,
                                       comb.fixed = FALSE)

forest(intensity.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(habitat.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(firetype.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(region.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(biome.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)

# Funnel Plots

#abundance
funnel(m.abundance.hksj, 
       xlab="Hedges' g",
       contour = c(0.95, 0.975, 0.99),
       col.contour = c("grey25", "grey31", "grey61"))
legend(0.25, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
       fill=c("grey25", "grey31", "grey61"))

#richness
funnel(m.richness.hksj, 
       xlab="Hedges' g",
       contour = c(0.95, 0.975, 0.99),
       col.contour = c("grey25", "grey31", "grey61"))
legend(0.25, 0, c("p < 0.05", "p < 0.025", "p < 0.01"),bty = "n",
       fill=c("grey25", "grey31", "grey61"))

##############
## CARABIDS ##
##############

# pooled effect size using Hartung-Knapp-Sidik-Jonkman Method #
# abundance
m.abundance.hksj <- metagen(EffectSize, # effect sizes
                  EffectSizeStandardError, # effect size standard errors
                  data = carabidData[carabidData$Metric=="Abundance",], # data
                  studlab = paste(PaperID), # paper ID for random effect
                  comb.fixed = FALSE, # not doing fixed effects
                  comb.random = TRUE, # doing random effects
                  method.tau = "ML", # Sidik-Jonkman estimator for variance
                  hakn = TRUE, # HKSJ Method
                  prediction = TRUE, # predict an interval
                  sm="SMD") # mean difference effect size

# richness
m.richness.hksj <- metagen(EffectSize, # effect sizes
                            EffectSizeStandardError, # effect size standard errors
                            data = carabidData[carabidData$Metric=="Richness",], # data
                            studlab = paste(PaperID), # paper ID for random effect
                            comb.fixed = FALSE, # not doing fixed effects
                            comb.random = TRUE, # doing random effects
                            method.tau = "ML", # Sidik-Jonkman estimator for variance
                            hakn = TRUE, # HKSJ Method
                            prediction = TRUE, # predict an interval
                            sm="SMD") # mean difference effect size

# search for outliers and influence #
spot.outliers.random(data=m.abundance.hksj)
spot.outliers.random(data=m.richness.hksj)

m.abundance.hksj <- update.meta(m.abundance.hksj,
                                subset = PaperID %!in% c("Castillo and Wagner (2 years post)",
                                                         "Castillo and Wagner (4 years post)"))

m.richness.hksj <- update.meta(m.richness.hksj,
                               subset = PaperID %!in% c("Castillo and Wagner (2 years post)",
                                                        "Castillo and Wagner (4 years post)",
                                                        "Castillo and Wagner (3 years post)"))

# forest plots #

forest(m.abundance.hksj, digits.sd = 2, digits.se = 2)
forest(m.richness.hksj, digits.sd = 2, digits.se = 2)

# subgroup analysis

# abundance
intensity.abundance.subgroup <- update.meta(m.abundance.hksj,
                                  byvar = Intensity,
                                  comb.random = TRUE,
                                  comb.fixed = FALSE)
habitat.abundance.subgroup <- update.meta(m.abundance.hksj,
                               byvar = HabitatType,
                               comb.random = TRUE,
                               comb.fixed = FALSE)
firetype.abundance.subgroup <- update.meta(m.abundance.hksj,
                                          byvar = FireType,
                                          comb.random = TRUE,
                                          comb.fixed = FALSE)
region.abundance.subgroup <- update.meta(m.abundance.hksj,
                                         byvar = Region,
                                         comb.random = TRUE,
                                         comb.fixed = FALSE)
biome.abundance.subgroup <- update.meta(m.abundance.hksj,
                                         byvar = Biome,
                                         comb.random = TRUE,
                                         comb.fixed = FALSE)

forest(intensity.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(habitat.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(firetype.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(region.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(biome.abundance.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)

# richness
intensity.richness.subgroup <- update.meta(m.richness.hksj,
                                            byvar = Intensity,
                                            comb.random = TRUE,
                                            comb.fixed = FALSE)
habitat.richness.subgroup <- update.meta(m.richness.hksj,
                                          byvar = HabitatType,
                                          comb.random = TRUE,
                                          comb.fixed = FALSE)
firetype.richness.subgroup <- update.meta(m.richness.hksj,
                                           byvar = FireType,
                                           comb.random = TRUE,
                                           comb.fixed = FALSE)
region.richness.subgroup <- update.meta(m.richness.hksj,
                                         byvar = Region,
                                         comb.random = TRUE,
                                         comb.fixed = FALSE)
biome.richness.subgroup <- update.meta(m.richness.hksj,
                                        byvar = Biome,
                                        comb.random = TRUE,
                                        comb.fixed = FALSE)

forest(intensity.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(habitat.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(firetype.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(region.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)
forest(biome.richness.subgroup, layout = "subgroup", digits.sd = 2, digits.se = 2)

# Funnel Plots

#abundance
funnel(m.abundance.hksj, 
       xlab="Hedges' g",
       contour = c(0.95, 0.975, 0.99),
       col.contour = c("grey25", "grey31", "grey61"))
legend(0.25, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
       fill=c("grey25", "grey31", "grey61"))

#richness
funnel(m.richness.hksj, 
       xlab="Hedges' g",
       contour = c(0.95, 0.975, 0.99),
       col.contour = c("grey25", "grey31", "grey61"))
legend(0.25, 0, c("p < 0.05", "p < 0.025", "p < 0.01"),bty = "n",
       fill=c("grey25", "grey31", "grey61"))

#################
## BUTTERFLIES ##
#################

##############
## COMBINED ##
##############