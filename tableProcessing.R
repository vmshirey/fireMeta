#####################################################################################
## Script for calculating species richness and species abundance from tabular data ##
## Author: Vaughn M. Shirey                                                        ##
#####################################################################################

paper <- read.csv("Paper_001.csv", header=TRUE, sep=",")

head(paper)

#######################
## species abundance ##
#######################

sum(paper[ which(paper$Site==2014),]$HN) 

sum(paper[ which(paper$Site==2013),]$L)
sum(paper[ which(paper$Site==2014),]$L)

sum(paper[ which(paper$Site==2014),]$SB)

sum(paper[ which(paper$Site==2013),]$CCW)
sum(paper[ which(paper$Site==2015),]$CCW)

sum(paper[ which(paper$Site==2013),]$HW)
sum(paper[ which(paper$Site==2014),]$HW)

sum(paper[ which(paper$Site==2013),]$CCE)
sum(paper[ which(paper$Site==2014),]$CCE)

sum(paper[ which(paper$Site==2013),]$SF)
sum(paper[ which(paper$Site==2014),]$SF)
sum(paper[ which(paper$Site==2015),]$SF)

sum(paper[ which(paper$Site==2013),]$HLP)
sum(paper[ which(paper$Site==2014),]$HLP)

sum(paper[ which(paper$Site==2013),]$WH)
sum(paper[ which(paper$Site==2015),]$WH)

sum(paper[ which(paper$Site==2013),]$MU)
sum(paper[ which(paper$Site==2014),]$MU)

sum(paper[ which(paper$Site==2014),]$MR)

######################
## species richness ##
######################

sum(paper[ which(paper$Site==2014),]$HN != 0) 

sum(paper[ which(paper$Site==2013),]$L != 0)
sum(paper[ which(paper$Site==2014),]$L != 0)

sum(paper[ which(paper$Site==2014),]$SB != 0)

sum(paper[ which(paper$Site==2013),]$CCW != 0)
sum(paper[ which(paper$Site==2015),]$CCW != 0)

sum(paper[ which(paper$Site==2013),]$HW != 0)
sum(paper[ which(paper$Site==2014),]$HW != 0)

sum(paper[ which(paper$Site==2013),]$CCE != 0)
sum(paper[ which(paper$Site==2014),]$CCE != 0)

sum(paper[ which(paper$Site==2013),]$SF != 0)
sum(paper[ which(paper$Site==2014),]$SF != 0)
sum(paper[ which(paper$Site==2015),]$SF != 0)

sum(paper[ which(paper$Site==2013),]$HLP != 0)
sum(paper[ which(paper$Site==2014),]$HLP != 0)

sum(paper[ which(paper$Site==2013),]$WH != 0)
sum(paper[ which(paper$Site==2014),]$WH != 0)

sum(paper[ which(paper$Site==2013),]$MU != 0)
sum(paper[ which(paper$Site==2014),]$MU != 0)

sum(paper[ which(paper$Site==2014),]$MR != 0)