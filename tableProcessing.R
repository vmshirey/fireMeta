#####################################################################################
## Script for calculating species richness and species abundance from tabular data ##
## Author: Vaughn M. Shirey                                                        ##
#####################################################################################

require(compute.es)

paper <- read.csv("Carabidae/Paper_001.csv", header=TRUE, sep=",")
head(paper)

## Sites with fire
## Controls 2013: AG, HN, SB, FC, TC
## Controls 2014: AG, HF, FC, TC
## Controls 2015: AG, FC, TC

## Y02013: L, CCW, HW, CCE, SF, HLP, WH, MU
## Y02014: HN, SB

## Y12014: CCW, WH
## Y12015: SB, MR

###########################################################################################################################################
##
## Control2013: AG, HN, SB, FC, TC
## Control2014: AG, HF, FC, TC
## Control2015: AG, FC, TC
##
## Y02013: L, CCW, HW, CCE, SF, HLP, WH, MU
## Y02014: HN, SB | Repeat: L, HW, CCE, SF, HLP, MU
## Y02015: CCW, WH | Repeats: SF ---- Trash
##
## Y12014: CCW, WH
## Y12015: SB, MR | Repeat: CCE, MU
##
###########################################################################################################################################

################
## Abundances ##
################

controls2013 <- c(sum(paper[ which(paper$Site==2013),]$AG), sum(paper[ which(paper$Site==2013),]$HN), sum(paper[ which(paper$Site==2013),]$SB),
                  sum(paper[ which(paper$Site==2013),]$FC), sum(paper[ which(paper$Site==2013),]$TC))
controls2014 <- c(sum(paper[ which(paper$Site==2014),]$AG), sum(paper[ which(paper$Site==2014),]$HF), sum(paper[ which(paper$Site==2014),]$FC),
                  sum(paper[ which(paper$Site==2014),]$TC))
controls2015 <- c(sum(paper[ which(paper$Site==2015),]$AG), sum(paper[ which(paper$Site==2015),]$FC),
                  sum(paper[ which(paper$Site==2015),]$TC))

print(paste("Control 2013 Mean: ", mean(controls2013), " | Control 2013 SE: ", sd(controls2013)/sqrt(length(controls2013))))
print(paste("Control 2014 Mean: ", mean(controls2014), " | Control 2014 SE: ", sd(controls2013)/sqrt(length(controls2014))))
print(paste("Control 2015 Mean: ", mean(controls2015), " | Control 2015 SE: ", sd(controls2013)/sqrt(length(controls2015))))

y02013 <- c(sum(paper[ which(paper$Site==2013),]$L), sum(paper[ which(paper$Site==2013),]$CCW), sum(paper[ which(paper$Site==2013),]$HW),
                  sum(paper[ which(paper$Site==2013),]$CCE), sum(paper[ which(paper$Site==2013),]$SF), sum(paper[ which(paper$Site==2013),]$HLP),
                  sum(paper[ which(paper$Site==2013),]$WH), sum(paper[ which(paper$Site==2013),]$MU))
y02014 <- c(sum(paper[ which(paper$Site==2014),]$SB), sum(paper[ which(paper$Site==2014),]$HN))
re2014 <- c(sum(paper[ which(paper$Site==2014),]$L),
                  sum(paper[ which(paper$Site==2014),]$HW), sum(paper[ which(paper$Site==2014),]$CCE), sum(paper[ which(paper$Site==2014),]$SF),
                  sum(paper[ which(paper$Site==2014),]$HLP), sum(paper[ which(paper$Site==2014),]$MU))
y02015 <- c(sum(paper[ which(paper$Site==2015),]$CCW), sum(paper[ which(paper$Site==2015),]$WH))
re2015 <- sum(paper[ which(paper$Site==2015),]$SF)

print(paste("Y0 2013 Mean: ", mean(y02013), " | Y0 2013 SE: ", sd(y02013)/sqrt(length(y02013))))
print(paste("Y0 2014 Mean: ", mean(y02014), " | Y0 2014 SE: ", sd(y02014)/sqrt(length(y02014))))
print(paste("RE 2014 Mean: ", mean(re2014), " | RE 2014 SE: ", sd(re2014)/sqrt(length(re2014))))
print(paste("Y0 2015 Mean: ", mean(y02015), " | Y0 2015 SE: ", sd(y02015)/sqrt(length(y02015))))
print(paste("RE 2015 Mean: ", mean(re2015), " | RE 2015 SE: ", sd(re2015)/sqrt(length(re2015))))

y12014 <- c(sum(paper[ which(paper$Site==2014),]$CCW), sum(paper[ which(paper$Site==2014),]$WH))
y12015 <- c(sum(paper[ which(paper$Site==2015),]$SB), sum(paper[ which(paper$Site==2015),]$MR))
re12015 <- c(sum(paper[ which(paper$Site==2015),]$CCE), sum(paper[ which(paper$Site==2015),]$MU))

print(paste("Y1 2014 Mean: ", mean(y12014), " | Y1 2014 SE: ", sd(y12014)/sqrt(length(y12014))))
print(paste("Y1 2015 Mean: ", mean(y12015), " | Y1 2015 SE: ", sd(y12015)/sqrt(length(y12015))))
print(paste("RE1 2015 Mean: ", mean(re12015), " | RE1 2015 SE: ", sd(re12015)/sqrt(length(re12015))))

###############
## Hedge's d ##
###############

mes(mean(controls2013), mean(y02013), sd(controls2013), sd(y02013), length(controls2013), length(y02013))
mes(mean(controls2014), mean(y02014), sd(controls2014), sd(y02014), length(controls2014), length(y02014))
mes(mean(controls2014), mean(re2014), sd(controls2014), sd(re2014), length(controls2014), length(re2014))
mes(mean(controls2015), mean(y02015), sd(controls2015), sd(y02015), length(controls2015), length(y02015))
mes(mean(controls2015), mean(re2015), sd(controls2015), sd(re2015), length(controls2015), length(re2015))

mes(mean(controls2014), mean(y12014), sd(controls2014), sd(y12014), length(controls2014), length(y12014))
mes(mean(controls2015), mean(y12015), sd(controls2015), sd(y12015), length(controls2015), length(y12015))
mes(mean(controls2015), mean(re12015), sd(controls2015), sd(re12015), length(controls2015), length(re12015))

###########################################################################################################################################
##
## Control2013: AG, HN, SB, FC, TC
## Control2014: AG, HF, FC, TC
## Control2015: AG, FC, TC
##
## Y02013: L, CCW, HW, CCE, SF, HLP, WH, MU
## Y02014: HN, SB | Repeat: L, HW, CCE, SF, HLP, MU
## Y02015: CCW, WH | Repeats: SF
##
## Y12014: CCW, WH
## Y12015: SB, MR | Repeat: CCE, MU
##
###########################################################################################################################################

controls2013 <- c(sum(paper[ which(paper$Site==2013),]$AG !=0 ), sum(paper[ which(paper$Site==2013),]$HN != 0), sum(paper[ which(paper$Site==2013),]$SB != 0),
                  sum(paper[ which(paper$Site==2013),]$FC != 0), sum(paper[ which(paper$Site==2013),]$TC != 0))
controls2014 <- c(sum(paper[ which(paper$Site==2014),]$AG != 0), sum(paper[ which(paper$Site==2014),]$HF != 0), sum(paper[ which(paper$Site==2014),]$FC != 0),
                  sum(paper[ which(paper$Site==2014),]$TC != 0))
controls2015 <- c(sum(paper[ which(paper$Site==2015),]$AG != 0), sum(paper[ which(paper$Site==2015),]$FC != 0),
                  sum(paper[ which(paper$Site==2015),]$TC != 0))

print(paste("Control 2013 Mean: ", mean(controls2013), " | Control 2013 SE: ", sd(controls2013)/sqrt(length(controls2013))))
print(paste("Control 2014 Mean: ", mean(controls2014), " | Control 2014 SE: ", sd(controls2013)/sqrt(length(controls2014))))
print(paste("Control 2015 Mean: ", mean(controls2015), " | Control 2015 SE: ", sd(controls2013)/sqrt(length(controls2015))))

y02013 <- c(sum(paper[ which(paper$Site==2013),]$L != 0), sum(paper[ which(paper$Site==2013),]$CCW != 0), sum(paper[ which(paper$Site==2013),]$HW != 0),
            sum(paper[ which(paper$Site==2013),]$CCE != 0), sum(paper[ which(paper$Site==2013),]$SF != 0), sum(paper[ which(paper$Site==2013),]$HLP != 0),
            sum(paper[ which(paper$Site==2013),]$WH != 0), sum(paper[ which(paper$Site==2013),]$MU != 0))
y02014 <- c(sum(paper[ which(paper$Site==2014),]$SB != 0), sum(paper[ which(paper$Site==2014),]$HN != 0))
re2014 <- c(sum(paper[ which(paper$Site==2014),]$L != 0),
            sum(paper[ which(paper$Site==2014),]$HW != 0), sum(paper[ which(paper$Site==2014),]$CCE != 0), sum(paper[ which(paper$Site==2014),]$SF != 0),
            sum(paper[ which(paper$Site==2014),]$HLP != 0), sum(paper[ which(paper$Site==2014),]$MU != 0))
y02015 <- c(sum(paper[ which(paper$Site==2015),]$CCW != 0), sum(paper[ which(paper$Site==2015),]$WH != 0))
re2015 <- sum(paper[ which(paper$Site==2015),]$SF != 0)

print(paste("Y0 2013 Mean: ", mean(y02013), " | Y0 2013 SE: ", sd(y02013)/sqrt(length(y02013))))
print(paste("Y0 2014 Mean: ", mean(y02014), " | Y0 2014 SE: ", sd(y02014)/sqrt(length(y02014))))
print(paste("RE 2014 Mean: ", mean(re2014), " | RE 2014 SE: ", sd(re2014)/sqrt(length(re2014))))
print(paste("Y0 2015 Mean: ", mean(y02015), " | Y0 2015 SE: ", sd(y02015)/sqrt(length(y02015))))
print(paste("RE 2015 Mean: ", mean(re2015), " | RE 2015 SE: ", sd(re2015)/sqrt(length(re2015))))

y12014 <- c(sum(paper[ which(paper$Site==2014),]$CCW != 0), sum(paper[ which(paper$Site==2014),]$WH != 0))
y12015 <- c(sum(paper[ which(paper$Site==2015),]$SB != 0), sum(paper[ which(paper$Site==2015),]$MR != 0))
re12015 <- c(sum(paper[ which(paper$Site==2015),]$CCE != 0), sum(paper[ which(paper$Site==2015),]$MU != 0))

print(paste("Y1 2014 Mean: ", mean(y12014), " | Y1 2014 SE: ", sd(y12014)/sqrt(length(y12014))))
print(paste("Y1 2015 Mean: ", mean(y12015), " | Y1 2015 SE: ", sd(y12015)/sqrt(length(y12015))))
print(paste("RE1 2015 Mean: ", mean(re12015), " | RE1 2015 SE: ", sd(re12015)/sqrt(length(re12015))))

###############
## Hedge's d ##
###############

mes(mean(controls2013), mean(y02013), sd(controls2013), sd(y02013), length(controls2013), length(y02013))
mes(mean(controls2014), mean(y02014), sd(controls2014), sd(y02014), length(controls2014), length(y02014))
mes(mean(controls2014), mean(re2014), sd(controls2014), sd(re2014), length(controls2014), length(re2014))
mes(mean(controls2015), mean(y02015), sd(controls2015), sd(y02015), length(controls2015), length(y02015))
mes(mean(controls2015), mean(re2015), sd(controls2015), sd(re2015), length(controls2015), length(re2015))

mes(mean(controls2014), mean(y12014), sd(controls2014), sd(y12014), length(controls2014), length(y12014))
mes(mean(controls2015), mean(y12015), sd(controls2015), sd(y12015), length(controls2015), length(y12015))
mes(mean(controls2015), mean(re12015), sd(controls2015), sd(re12015), length(controls2015), length(re12015))
