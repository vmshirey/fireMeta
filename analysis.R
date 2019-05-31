## Analysis code for pyroentomology meta-analysis
## Vaughn Shirey, 2019

# load libraries
library(meta)
library(metafor)

# read in data
carabidData <- read.csv("carabidData.csv", header=TRUE)

## BEES ##

## CARABIDS ##

# pooled effect size using Hartung-Knapp-Sidik-Jonkman Method
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

# forest plots
forest(m.abundance.hksj, digits.sd = 2, digits.se = 2)
forest(m.richness.hksj, digits.sd = 2, digits.se = 2)


# subgroup analysis
intensity.abundance.subgroup <- update.meta(m.abundance.hksj,
                                  byvar = Intensity,
                                  comb.random = TRUE,
                                  comb.fixed = FALSE)

forest(intensity.abundance.subgroup)

habitat.abundance.subgroup <- update.meta(m.abundance.hksj,
                               byvar = HabitatType,
                               comb.random = TRUE,
                               comb.fixed = FALSE)

forest(habitat.abundance.subgroup)

firetype.abundance.subgroup <- update.meta(m.abundance.hksj,
                                          byvar = FireType,
                                          comb.random = TRUE,
                                          comb.fixed = FALSE)

forest(firetype.abundance.subgroup)

region.abundance.subgroup <- update.meta(m.abundance.hksj,
                                         byvar = Region,
                                         comb.random = TRUE,
                                         comb.fixed = FALSE)

forest(region.abundance.subgroup)

intensity.richness.subgroup <- update.meta(m.richness.hksj,
                                            byvar = Intensity,
                                            comb.random = TRUE,
                                            comb.fixed = FALSE)

forest(intensity.richness.subgroup)

habitat.richness.subgroup <- update.meta(m.richness.hksj,
                                          byvar = HabitatType,
                                          comb.random = TRUE,
                                          comb.fixed = FALSE)

forest(habitat.richness.subgroup)

firetype.richness.subgroup <- update.meta(m.richness.hksj,
                                           byvar = FireType,
                                           comb.random = TRUE,
                                           comb.fixed = FALSE)

forest(firetype.richness.subgroup)

region.richness.subgroup <- update.meta(m.richness.hksj,
                                         byvar = Region,
                                         comb.random = TRUE,
                                         comb.fixed = FALSE)

forest(region.richness.subgroup)

## BUTTERFLIES ##