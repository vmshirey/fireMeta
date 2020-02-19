#####################################################################################
## Script for calculating effect sizes from different statistical tests.           ##
## Author: Vaughn M. Shirey                                                        ##
#####################################################################################

## 95% CI to SD = sqrt(n)*(upper - lower)/3.92

library(compute.es)
library(esc)

################
## Vogel 2007 ##
################
mes(31.48, 27.78, 7.25, 5.67, ) # abundance
mes() # richness

###############
## Riva 2019 ##
###############
mes(16.06, 15.5, 5.3, 5.6, 18, 12) # richness 
mes(73.33, 86.58, 49.33, 52.33, 18, 12) # abundance

#####################
## Rudolph and Ely ##
#####################
ub <- c(242.7, 57.6, 86.8, 153.0)
b <- c(348.1, 200.3, 260.5, 502.9)
mes(mean(b), mean(ub), sd(b), sd(ub), 4, 4) # abundance loblolly

ub <- c(110.0, 106.7, 45.8, 50.0)
b <- c(275.0, 165.0, 56.7, 126.7)
mes(mean(b), mean(ub), sd(b), sd(ub), 4, 4) # abundance loblolly

#######################
## Salvato & Salvato ##
#######################
esc_t(t=1.7885, p=0.0669, grp1n=29, grp2n=4, es.type="g") # abundance

ub <- c(31, 15, 8, 2, 6, 8)
b <- c(23, 14, 9, 2, 2, 8)
mes(mean(b), mean(ub), sd(b), sd(ub), 6, 6) # abundance

#################
## Ohwaki 2019 ##
#################
mes(10.91, 17.75, 3.62, 2.36, 11, 4) # richness
mes(25.27, 68.75, 13.09, 16.15, 11, 4) # abundance

############################
## Pascale and Theit 2006 ##
############################


############################
## McCullough et al. 2019 ##
############################
# butterflies #
mes(0.004384, 0.005591, 0.003882, 0.008501, 8, 12) # mid term burn abundance 2012
mes(0.000954, 0.005591, 0.001102, 0.008501, 4, 12) # low term burn abundance 2012

mes(0.003647, 0.00048, 0.006195, 0.00107, 57, 18) # mid term burn abundance 2014
mes(0.000989, 0.00048, 0.000857, 0.00107, 3, 18) # low term burn abundance 2014

mes(0.002039, 0.001485, 0.003514, 0.00217, 132, 72) # mid term burn abundance 2015
mes(0.000502, 0.001485, 0.001182, 0.00217, 24, 72) # low term burn abundance 2015

mes(0.0000924, 0.000456, 0.000316, 0.000877, 24, 30) # low term burn abundance 2016
mes(0.0000378, 0.000456, 0.00128, 0.000877, 84, 30) # mid term burn abundance 2016


##################
## Gaigher 2019 ##
##################
# thicket #
mes(71.5, 95.75, 56.71, 60.15, 8, 4) # abundance
mes(6.88, 11.25, 1.25, 1.89, 8, 4) # richness

# 1-2 y #
mes(46, 95.75, 14.66, 60.15, 7, 4) # abundance
mes(12.57, 11.25, 1.13, 1.89, 7, 4) # richness

# < 1y #
mes(51.88, 95.75, 19.37, 60.15, 17, 4) # abundance
mes(12.88, 2.80, 11.25, 1.89, 17, 4) # richness

############
## Oliver ##
############
# ants #
mes() # richness

#################
## Kobach 2013 ##
#################
mes(0.625, 2.125, 1.19, 1.96, 7, 8) # abundance
mes(0.625, 1.88, 1.19, 1.73, 7, 8) # richness

##################
## Kambach 2013 ##
##################
# bees #
esc_t(t=-3.80, p=0.002, totaln=15, es.type="g") # richness
esc_t(t=6.05, p=0.001, totaln=15, es.type="g") # abundance

#################
## Doabma 2014 ##
#################
# carabidae #
mes()

##################################
## Valley Forge NHP Bee Surveys ##
##################################
# bees #
mes(42, 41, 0.00, 23.14, 1, 6) # 2014 burn abundance

########################
## Kwilosz  & Knutson ##
########################
# Butterflies #
control <- c(104, 65)
burn <- c(288, 213)
mes(mean(burn), mean(control), sd(burn), sd(control), 2, 2)

# Freq 1995 #
control <- c(68, 131, 159, 121, 59, 109, 79, 35, 58, 37, 20, 26, 21, 20, 45, 26)
burn <- c(245, 323, 207, 126, 165, 215, 276)
mes(mean(burn), mean(control), sd(burn), sd(control), length(burn), length(control))

# Freq 1996 #
control <- c(6, 20, 15, 2, 44, 63, 89, 86, 5, 6, 13, 4, 0, 12, 23, 25, 17, 60, 33, 39, 25, 0, 5, 14, 88, 103)
burn <- c(52, 144, 39, 158, 224, 284, 335, 21, 30, 41, 29, 5, 16, 92, 161, 164)
mes(mean(burn), mean(control), sd(burn), sd(control), length(burn), length(control))

# Freq 1997 #
control <- c(60, 79, 36, 103, 96, 50, 26, 31, 23, 0, 71, 46, 37, 75, 71, 78, 164, 184, 107)
burn <- c(160, 153, 152, 365, 725, 541, 70, 94, 93, 24, 189, 296, 248)
mes(mean(burn), mean(control), sd(burn), sd(control), length(burn), length(control))

##################
## Moylett 2019 ##
##################
# Bees #
mes(5.97, 2.5, 2*(8.36-4.25)/3.92, 2*(3.88-1.72)/3.92, 4, 4) #0y abun
mes(4.14, 2.5, 2*(5.86-2.95)/3.92, 2*(3.88-1.72)/3.92, 4, 4) #1y abun
mes(4.10, 2.5, 2*(5.75-2.91)/3.92, 2*(3.88-1.72)/3.92, 4, 4) #2y abun

mes(4.02, 1.85, 2*(5.05-3.2)/3.92, 2*(2.46-1.44)/3.92, 4, 4) #0y rich
mes(3.27, 1.85, 2*(4.08-2.59)/3.92, 2*(2.46-1.44)/3.92, 4, 4) #1y rich
mes(3.02, 1.85, 2*(3.83-2.41)/3.92, 2*(2.46-1.44)/3.92, 4, 4) #2y rich

############
## Rubene ##
############
# Bees #
mes(19.2, 16.824, 6.03, 4.32, 15, 17) # richness
mes(250.133, 109.71, 179.65, 58.196, 15, 17) # abundance

##########################
## DeSantis et al. 2007 ##
##########################
# carabidae #
# fire #
mes(2.35, 5.05, 0.806225775, 1.857417562, 4, 4) # abun

# fire + mechanical #
mes(3.2, 5.05, 0.909212113, 1.857417562, 4, 4) # abun

############################
## Van Amburg et al. 1981 ##
############################
esc_chisq(chisq=8.33, p=0.01, totaln=240, es.type="g")

##################
## Adedoja 2019 ##
##################
mes(38, 24.83, 22.52, 15.33, 6, 3) # abundance
mes(10, 7.83, 3.61, 4.54, 6, 3) # richness

##########################
## McFarren et al. 1995 ##
##########################
# Carabids #
esc_f(f=7.92, grp1n = 19, grp2n = 19, es.type="g")

################
## Sasal 2015 ##
################
# Carabidae #
# Fire & Grazing, Control was Grazed #


#################
## Decker 2019 ##
#################
# Bees #
mes(76.86354, 39.63636, 88.79849, 33.00872, 22, 22) # y0 growing abundance
mes(13.18182, 10.86364, 9.912168, 7.845505, 22, 22) # y0 growing richness
mes(74.2381, 39.63636, 92.45048, 33.00872, 22, 21) # y0 dormant abundance
mes(12.66667, 10.86364, 9.520154, 7.845505, 22, 21) # y0 dormant richness

#################
## Wikers 1995 ##
#################
mes(46.8, 29.5, 9.8, 5.2, 5, 5) # burn abundance
mes(51.0, 29.5, 27.5, 5.2, 5, 5) # extra fuel burn abundance

#########################
## Villa-Castillo 2002 ##
#########################
## Fire ##
mes(9, 1.75, 1.83, 0.95, 4, 4) # 1998 richness
mes(8, 2.25, 2.44, 1.89, 4, 4) # 1999 richness
mes(8.25, 1.25, 2.98, 0.43, 4, 4) # 2000 richness

mes(145, 122, 50.61, 89.01, 4, 4) # 1998 abundance
mes(49.25, 151, 22.44, 101.96, 4, 4) # 1999 abundance
mes(25, 13.5, 7.39, 3.77, 4, 4) # 2000 abundance

## Mechanical + Fire ##
mes(5.75, 1.75, 0.95, 0.95, 4, 4) # 1998 richness
mes(5, 2.25, 1.41, 1.89, 4, 4) # 1999 richness
mes(3.75, 1.25, 1.70, 0.43, 4, 4) # 2000 richness

mes(55.25, 122, 5.85, 89.01, 4, 4) # 1998 abundance
mes(23.25, 151, 16.54, 101.96, 4, 4) # 1999 abundance
mes(25.50, 13.5, 23.04, 3.77, 4, 4) # 2000 abundance

############
## Harris ##
############
mes(9.09, 44.45, 4.8, 19.07, 13, 13) #y0 abundance
mes(1.85, 5.62, 0.8, 1.15, 13, 13) #y0 richness

mes(17.0, 77.28, 15.98, 25.02, 17, 17) #y1 abundance
mes(2.24, 6.47, 2.12, 1.62, 17, 17) #y1 richness

###################
## Sackmann 2006 ##
###################
mes(305, 270, 435, 301, 8, 8) # abundance forest
mes(29.6, 60, 24.9, 60.9, 6, 6) # abundance scrub
mes(32.7, 16.5, 28.2, 9.9, 5, 5) # abundance steppe

##################
## Rickard 1970 ##
##################
c.abun <- c(20.8, 9.5, 1.7, 1.7, 41.8, 18.8, 5.2, 4.4, 6.3, 7.8, 3.2, 3.2, 4.6, 21.1, 7.8, 4.6, 2.7)
b.abun <- c(48.4, 18.3, 19.2, 19.2, 30.8, 8.7, 3.3, 4.4, 0, 0, 0, 1.0, 7.3, 13.8, 1.9, 4.6, 0)
mes(mean(b.abun), mean(c.abun), sd(b.abun), sd(c.abun), 17, 17)

#######################
## Muona and Rutanen ##
#######################
c.abun <- c()
b.abun <- c()
mes() # abundance

##########################
## Holiday 1991 & 1992 ##
##########################
# aspen forest
mes(22.25, 63.19, 1.77*sqrt(16), 4.78*sqrt(16), 16, 16) # abundance

# conifer forest
mes(17.31, 51.75, 1.62*sqrt(16), 6.30*sqrt(16), 16, 16) # abundance

# pooled #
esc_f(f=0.1, grp1n = 16, grp2n = 16, es.type="g") # 1992 paper richness

############
## Garcia ##
############
u.abun <- c(1.0, 0.2)
b.abun <- c(0.1, 0.1, 0.7)

mes(mean(b.abun), mean(u.abun), sd(b.abun), sd(u.abun), 10, 10)

##########################
## Dissertation Stanton ##
##########################
b.abun <- c(22.86, 20.42)
b.abun.sd <- c((7.72*sqrt(12))^2, (1.22*sqrt(12))^2)
mes(mean(b.abun), 18.85, sum(b.abun.sd)/12, 5.42, 12, 12) # no freq 0y

mes(21.99, 20.24, 1.56, 2.10, 12, 12) # no freq 1y

mes(20.24, 18.59, 3.32, 2.1, 12, 12) # no freq 2y

b.abun <- c(9.09, 7.61)
b.abun.sd <- c((1.58*sqrt(12))^2, (1.56*sqrt(12))^2)
mes(mean(b.abun), 9.88, sum(b.abun.sd)/12, 2.26, 12, 12) # freq

# abundance #
b.abun <- c(484.4, 554.4)
b.abun.sd <- c((196*sqrt(12))^2, (106.4*sqrt(12))^2)
mes(mean(b.abun), 403.2, sum(b.abun.sd)/12, 117.6, 12, 12) # no freq 0y

mes(280, 266, 89.6, 106.4, 12, 12) # no freq 1y

mes(184.8, 176.4, 55.2, 78.4, 12, 12) # no freq 2y

b.abun <- c(39.2, 39.2)
b.abun.sd <- c((16.8*sqrt(12))^2, (16.8*sqrt(12))^2)
mes(mean(b.abun), 30.8, sum(b.abun.sd)/12, 23.2, 12, 12) # freq

###################
## Hanula et al. ##
###################
c.abun <- c(44.17, 79.17)
b.abun <- c(88.2, 42.6)

c.abun.sd <- c((6.89*sqrt(6))^2, (8.44*sqrt(6))^2)
b.abun.sd <- c((8.3*sqrt(6))^2, (3.3*sqrt(6))^2)

mes(mean(b.abun), mean(c.abun), sum(b.abun.sd)/12, sum(c.abun.sd)/12, 12, 12)

##########################
## Campbell et. al 2018 ##
##########################
c.abun <- c(0.67, 5.67)
b.abun <- c(3.67, 3.00)

c.abun.sd <- c((0.33*sqrt(11))^2, (2.03*sqrt(11))^2)
b.abun.sd <- c((3.67*sqrt(11))^2, (3.00*sqrt(11))^2)

mes(mean(b.abun), mean(c.abun), sum(b.abun.sd)/2, sum(c.abun.sd)/2, 11, 11)

# Bees #
c.abun <- c(0.8, 8.0)
mb.abun <- c(7.6, 8.4)

c.abun.sd <- c((0.2*sqrt(11))^2, (2.5*sqrt(11))^2)
mb.abun.sd <- c((1.9*sqrt(11))^2, (1.6*sqrt(11))^2)

mes(mean(mb.abun), mean(c.abun), sum(mb.abun.sd)/2, sum(c.abun.sd)/2, 11, 11)

# Butterflies #
# Fire #
c.abun <- c(0.1, 0.04)
mb.abun <- c(0.6, 0.2)
b.abun <- c(0.3, 0.2)

c.abun.sd <- c((0.03*sqrt(11))^2, (0.02*sqrt(11))^2)
mb.abun.sd <- c((0.2*sqrt(11))^2, (0.1*sqrt(11))^2)
b.abun.sd <- c((0.09*sqrt(11))^2, (0.09*sqrt(11))^2)

mes(mean(b.abun), mean(c.abun), sum(b.abun.sd)/2, sum(c.abun.sd)/2, 11, 11)
mes(mean(mb.abun), mean(c.abun), sum(mb.abun.sd)/2, sum(c.abun.sd)/2, 11, 11)

#################
## Rice et al. ##
#################
b.abun <- c(0, 0, 0, 40)
c.abun <- c(40, 40, 20, 20)

mes(mean(b.abun), mean(c.abun), sd(b.abun), sd(c.abun), 4, 4) # abundance

###################
## Barber et al. ##
###################
# Carabidae #
# fire
mes(57.67, 70.71, 12.1, 64.91, 3, 7) # abun
mes(12.67, 9, 4.93, 4.83, 3, 7) # rich

# fire + grazing
mes(143.75, 70.71, 151.27, 64.91, 4, 7) # abun
mes(10.5, 9, 6.24, 4.83, 4, 7) # rich

##################
## Valko et. al ##
##################

mes(3.53, 3.10, 2.90, 2.90, 60, 60) # abundance
mes(2.03, 1.80, 1.25, 1.10, 60, 60) # richness

#################
## Bess et. al ##
#################

#######################################################################################
## Data was lifted from figure 3, on page 780. Control is all pre-fire points; burn  ##
## is all post-fire points (including fire).                                         ##
## Used DataThief, rounded to nearest integer value for richness and hundreths-place ##
## for abundance.                                                                    ##
#######################################################################################

control.abd <- c(8, 6, 13, 41, 676, 29, 8)
control.rch <- c(1, 3, 2, 1, 3, 1, 1)
  
wf.abd <-c(35, 4, 7, 71)
wf.rch <- c(3, 3, 3, 1)

mes(mean(wf.rch), mean(control.rch), sd(wf.rch), sd(control.rch), length(wf.rch), length(control.rch))
mes(mean(wf.abd), mean(control.abd), sd(wf.abd), sd(control.abd), length(wf.abd), length(control.abd))

######################
## Ghandi et al. ##
######################

######################
## Data in Figure 5 ##
######################

mes(7.5, 19, 8.5, 7, 12, 18) #abundance
mes(20, 16, 8, 8, 12, 12) #abundance

######################
## Fernandez et al. ##
######################

#####################
## Data in Table 1 ##
#####################

wf.abd <- c(193, 87, 18)
wf.rch <- c(10, 6, 3)
  
control.abd <- c(69, 20, 6)
control.rch <- c(8, 6, 3)

print(paste("Control Rch Mean: ", mean(control.rch), " | Control Rch SE: ", sd(control.rch)/sqrt(length(control.rch))))
print(paste("Control Abd Mean: ", mean(control.abd), " | Control Abd SE: ", sd(control.abd)/sqrt(length(control.abd))))

print(paste("WF Rch Mean: ", mean(wf.rch), " | WF Abd SE: ", sd(wf.rch)/sqrt(length(wf.rch))))
print(paste("WF Abd Mean: ", mean(wf.abd), " | WF Abd SE: ", sd(wf.abd)/sqrt(length(wf.abd))))

# Submitting to the MES function:
mes(mean(wf.rch), # average fire richness
    mean(control.rch), # average control richness
    sd(wf.rch), # standard deviation of fire richness
    sd(control.rch), # standard deviation of control richness
    length(wf.rch), # number of samples in fire sites
    length(control.rch)) # number of samples control sites

mes(mean(wf.abd), mean(control.abd), sd(wf.abd), sd(control.abd), length(wf.abd), length(control.abd))

#########################
## Moretti et al. 2004 ##
#########################
## Carabidae in Figure 2 ##
mes(8, 8, 0.55*sqrt(24), 0.5*sqrt(24), 24, 24) # richness
mes(90, 95, 30*sqrt(24), 15*sqrt(24), 24, 24) # abundance

## Bees in ANOVA, Figure 1 (only fire freq) ##
esc_f(f=7.23, grp1n=24, grp2n=24, es.type="g")

####################
## Nunes et al.   ##
####################

###################
## Data in Fig 1 ##
###################

mes(0.85, 0.8, 0.17*sqrt(700), 0.3*sqrt(700), 700, 700)
mes(0.92, 0.81, 0.5*sqrt(350), 0.6*sqrt(350), 350, 350)

mes(0.045, 0.047, 0.005*sqrt(700), 0.0045*sqrt(700), 700, 700)
mes(0.06, 0.08, 0.06*sqrt(350), 0.09*sqrt(350), 350, 350)

mes(0.2, 0.5, 0.1*sqrt(700), 0.2*sqrt(700), 700, 700)
mes(0.5, 0.5, 0.3*sqrt(350), 0.1*sqrt(350), 350, 350)

mes(0.045, 0.075, 0.03*sqrt(700), 0.025*sqrt(700), 700, 700)
mes(0.01, 0.015, 0.01*sqrt(350), 0.01*sqrt(350), 350, 350)

############################
## Saint-Germain et al.   ##
############################
#####################
## Data in Table 1 ##
#####################
controlsN <- c(0.14, 0.15, 13.78, 4.6, 0.06, 0.06, 0.06, .69, 6.82, 1.15, 0.20, 0.12, 2.25, 1.86, 0.11)
controlsSD <- c(0.41, 0.65, 56.97, 12.58, 0.24, 0.24, 0.24, 1.36, 12.37, 1.55, 0.85, 0.36, 7.15, 6.33, 0.47)
burnedN <- c(0.08, 0.09, 0.08, 0.66, 0.21, 0.08, 0.25, 0.13, 0.17, 0.08, 0.09, 3.37, 0.65, 0.60, 0.09, 0.09, 0.38)
burnedSD <- c(0.29, 0.31, 0.29, 1.34, 0.50, 0.29, 0.45, 0.43, 0.40, 0.29, 0.31, 3.04, 0.98, 1.30, 0.31, 0.31, 0.90)
controlsM <- cbind(controlsN, controlsSD)
burnedM <- cbind(burnedN, burnedSD)
controlsPSD <- sqrt(
  ((controlsM[1,1] - 1)*controlsM[1,2]^2 + 
  (controlsM[2,1] - 1)*controlsM[2,2]^2 +
  (controlsM[3,1] - 1)*controlsM[3,2]^2 +
  (controlsM[4,1] - 1)*controlsM[4,2]^2 +
  (controlsM[5,1] - 1)*controlsM[5,2]^2 +
  (controlsM[6,1] - 1)*controlsM[6,2]^2 +
  (controlsM[7,1] - 1)*controlsM[7,2]^2 +
  (controlsM[8,1] - 1)*controlsM[8,2]^2 +
  (controlsM[9,1] - 1)*controlsM[9,2]^2 +
  (controlsM[10,1] - 1)*controlsM[10,2]^2 +
  (controlsM[11,1] - 1)*controlsM[11,2]^2 +
  (controlsM[12,1] - 1)*controlsM[12,2]^2 +
  (controlsM[13,1] - 1)*controlsM[13,2]^2 +
  (controlsM[14,1] - 1)*controlsM[14,2]^2 +
  (controlsM[15,1] - 1)*controlsM[15,2]^2)/
    (sum(controlsM[,1])-nrow(controlsM))
)
burnedPSD <- sqrt(
  abs(((burnedM[1,1] - 1)*burnedM[1,2]^2 + 
     (burnedM[2,1] - 1)*burnedM[2,2]^2 +
     (burnedM[3,1] - 1)*burnedM[3,2]^2 +
     (burnedM[4,1] - 1)*burnedM[4,2]^2 +
     (burnedM[5,1] - 1)*burnedM[5,2]^2 +
     (burnedM[6,1] - 1)*burnedM[6,2]^2 +
     (burnedM[7,1] - 1)*burnedM[7,2]^2 +
     (burnedM[8,1] - 1)*burnedM[8,2]^2 +
     (burnedM[9,1] - 1)*burnedM[9,2]^2 +
     (burnedM[10,1] - 1)*burnedM[10,2]^2 +
     (burnedM[11,1] - 1)*burnedM[11,2]^2 +
     (burnedM[12,1] - 1)*burnedM[12,2]^2 +
     (burnedM[13,1] - 1)*burnedM[13,2]^2 +
     (burnedM[14,1] - 1)*burnedM[14,2]^2 +
     (burnedM[15,1] - 1)*burnedM[15,2]^2 +
     (burnedM[16,1] - 1)*burnedM[16,2]^2 +
     (burnedM[17,1] - 1)*burnedM[17,2]^2)/
        (sum(burnedM[,1])-nrow(burnedM)))
)
mes(sum(burnedM[,1]), sum(controlsM[,1]), burnedPSD, controlsPSD, 17, 15) # abundance
mes(11.58, 8.13, 13.89-11.58, 9.94-8.13, 95, 95) # richness

#################
## Samu et al. ##
#################
#####################
## Data in Table 2 ##
#####################
mes(103.6, 32, 38.89*sqrt(20), 14.58*sqrt(20), 20, 20) # abundance
mes(9.8, 8.4, 1.69*sqrt(20), 2.55*sqrt(20), 20, 20) # richness

########################
## Martikainen et al. ##
########################
# Carabids
# values are inverse effects (fire)
mes(103, 88, 19*sqrt(20), 34*sqrt(20), 20, 20) # abundance
mes(5.6, 5.38, 0.9*sqrt(20), (6.54-5.38)*sqrt(20), 20, 20) # richness

# fire + mechanical
mes(532.97, 103, 19*sqrt(20), (608.65-532.97)*sqrt(20), 20, 20) # 0m abun
mes(520.00, 103, 19*sqrt(20), (567.57-520.00)*sqrt(20), 20, 20) # 10m abun
mes(390.27, 103, 19*sqrt(20), (450.87-390.27)*sqrt(20), 20, 20) # 50m abun

mes(23.23, 5.38, (27.46-23.23)*sqrt(20), (6.54-5.38)*sqrt(20), 20, 20) # 0m rich
mes(25.77, 5.38, (29.15-25.77)*sqrt(20), (6.54-5.38)*sqrt(20), 20, 20) # 10m rich
mes(26.15, 5.38, (28.92-26.15)*sqrt(20), (6.54-5.38)*sqrt(20), 20, 20) # 50m rich

##########
## Niwa ##
##########
# Carabids
mes(553.6, 1090.1, 196.6*sqrt(9), 350.00*sqrt(9), 9, 9) # abundance
mes(8.4, 7.4, 0.6*sqrt(9), 0.5*sqrt(9), 9, 9) # richness

######################
## Parmenter et al. ##
######################

mes(7.63, 15.39, 4.995, 11.66, 32, 32) # abundance
mes(2.97, 1.45, 1.12, 0.50, 32, 32) # richness

#############
## Koivula ##
#############

mes(0.082, 0.0006, 0.085, 0.0017, 8, 8) # low intensity 1 year
mes(0.271, 0.0006, 0.226, 0.0017, 8, 8) # high intensity 1 year

mes(0.0237, 0.0012, 0.0329, 0.0034, 8, 8) # low intensity 2 year
mes(0.0522, 0.0012, 0.0836, 0.0034, 8, 8) # high intensity 2 year

mes(0.5, 0.125, 0.535, 0.354, 8, 8) # richness 2 year low
mes(0.625, 0.125, 0.518, 0.354, 8, 8) # richness 2 year high

###########
## Pryke ##
###########
cbb <- c(1,11,0)
cbu <- c(11, 15, 2, 3)
fb <- c(2,7,1)
fu <- c(1,2,5,8)

mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 4)
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

cbb <- c(1,5,0)
cbu <- c(3,4,1,1)
fb <- c(1,2,1)
fu <- c(1,2,4,1)

mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 4)
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

#####################
## Gongalsky 2006  ##
#####################
unburn <- c(0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            3,
            0,
            0,
            4,
            1,
            2,
            2,
            3,
            2,
            1,
            2,
            3,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0)
burn <- c(14,
          17,
          30,
          7,
          13,
          15,
          1,
          0,
          7,
          3,
          3,
          2,
          16,
          7,
          1,
          5,
          1,
          11)

mes(mean(burn), mean(unburn), sd(burn), sd(unburn), length(burn), length(unburn))

unburn <- c(0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            2,
            0,
            0,
            3,
            1,
            2,
            2,
            1,
            2,
            1,
            1,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0)
burn <- c(3,
          4,
          4,
          3,
          3,
          5,
          1,
          0,
          4,
          3,
          2,
          2,
          4,
          2,
          1,
          3,
          1,
          3)
mes(mean(burn), mean(unburn), sd(burn), sd(unburn), length(burn), length(unburn))

burn <- c(16,
          56,
          14,
          43,
          40,
          20,
          93,
          46,
          33,
          31,
          58,
          0,
          46,
          53,
          54,
          40,
          44,
          91,
          62,
          63,
          20,
          33,
          14,
          6,
          13,
          3,
          17,
          14,
          1,
          14,
          17,
          9,
          32,
          13,
          18,
          16,
          26,
          63,
          12,
          9,
          8,
          8,
          3,
          12,
          16,
          2,
          20,
          29,
          15,
          5,
          9,
          13,
          5,
          7,
          2,
          7,
          18,
          7,
          15,
          5,
          0,
          9,
          16,
          11,
          87,
          1,
          18,
          11,
          3,
          16,
          33,
          12,
          28,
          6,
          44,
          4,
          8,
          35,
          5,
          4)
unburn <- c(39,
            28,
            40,
            58,
            31,
            8,
            83,
            60,
            97,
            31,
            13,
            0,
            3,
            3,
            5,
            3,
            27,
            20,
            0,
            26,
            6,
            20,
            6,
            8,
            21,
            9,
            11,
            10,
            7,
            16,
            19,
            1,
            0,
            13,
            5,
            8,
            25,
            2,
            8,
            5,
            0,
            8,
            4,
            5,
            19,
            8,
            16,
            3,
            5,
            6,
            6,
            0,
            0,
            0,
            0,
            1,
            12,
            19,
            0,
            15,
            4,
            43,
            0,
            2,
            14,
            5,
            28,
            5,
            126,
            12)
mes(mean(burn), mean(unburn), sd(burn), sd(unburn), length(burn), length(unburn))

burn <- c(2,
          7,
          5,
          13,
          7,
          7,
          18,
          8,
          4,
          9,
          12,
          0,
          11,
          15,
          11,
          7,
          9,
          8,
          4,
          13,
          10,
          10,
          8,
          5,
          11,
          2,
          9,
          7,
          1,
          9,
          9,
          5,
          8,
          6,
          10,
          5,
          18,
          7,
          4,
          3,
          5,
          5,
          3,
          6,
          12,
          2,
          11,
          15,
          7,
          4,
          6,
          8,
          3,
          2,
          1,
          5,
          9,
          4,
          8,
          4,
          0,
          8,
          11,
          8,
          25,
          1,
          6,
          7,
          3,
          6,
          9,
          6,
          11,
          5,
          10,
          3,
          6,
          10,
          5,
          3)
unburn <- c(10,
            14,
            18,
            14,
            10,
            7,
            16,
            19,
            12,
            16,
            7,
            0,
            3,
            3,
            4,
            3,
            10,
            12,
            0,
            15,
            4,
            6,
            3,
            5,
            4,
            4,
            2,
            3,
            3,
            5,
            10,
            1,
            0,
            8,
            5,
            6,
            10,
            2,
            5,
            4,
            0,
            4,
            4,
            4,
            4,
            3,
            7,
            2,
            2,
            2,
            6,
            0,
            0,
            0,
            0,
            1,
            6,
            11,
            0,
            8,
            2,
            10,
            0,
            1,
            4,
            5,
            12,
            4,
            14,
            5)
mes(mean(burn), mean(unburn), sd(burn), sd(unburn), length(burn), length(unburn))

# fire + mechanical #
mes(23.3, 16.3, 30.424, 22.856, 70, 70) # clearcut abun
mes(6.843, 5.771, 4.128, 4.879, 70, 70) # clearcut rich

mes(34.063, 16.3, 26.648, 22.856, 70, 70) # altcut abun
mes(9.888, 5.771, 4.600, 4.879, 70, 70) # altcut rich

######################
## Gongalsky et al. ##
######################
y1.ub <- c(21)
y1.b <- c(100, 115, 163, 398, 105)
y2.ub <- c(25)
y2.b <- c(49, 37, 51, 145, 76)
mes(mean(y1.b), mean(y1.ub), sd(y1.b), sd(y1.ub), 1, 5) # y1 abun
mes(mean(y2.b), mean(y2.ub), sd(y2.b), sd(y2.ub), 1, 5) # y2 abun

y1.ub <- c(7)
y1.b <- c(9, 15, 13, 11, 14)
y2.ub <- c(12)
y2.b <- c(9, 6, 8, 11, 5)
mes(mean(y1.b), mean(y1.ub), sd(y1.b), sd(y1.ub), 1, 5) # y1 rich
mes(mean(y2.b), mean(y2.ub), sd(y2.b), sd(y2.ub), 1, 5) # y2 rich

#############
## Hammond ##
#############
burn <- c(24,
          14,
          17,
          11,
          16,
          25,
          2,
          6,
          4,
          23,
          21,
          13)

unburn <- c(51,
            81,
            122,
            115,
            129,
            38,
            25,
            51,
            59,
            74,
            60,
            71)

mes(mean(burn), mean(unburn), sd(burn), sd(unburn), length(burn), length(unburn))

burn <- c(6,
          4,
          5,
          6,
          6,
          7,
          2,
          4,
          3,
          3,
          7,
          4)
unburn <- c(7,
            10,
            8,
            8,
            8,
            9,
            5,
            6,
            7,
            7,
            7,
            6)
mes(mean(burn), mean(unburn), sd(burn), sd(unburn), length(burn), length(unburn))

#####################
## Toivenen et al. ##
#####################
# carabid #
# burn #
b.abun <- c(31, 62, 3, 4, 50, 102, 68, 41, 35)
u.abun <- c(28, 72, 89, 44, 100, 51, 28, 44, 43)

b.rich <- c(7, 10, 3, 3, 5, 10, 9, 7, 5)
u.rich <- c(5, 5, 8, 5, 8, 6, 9, 6, 8)

mes(mean(b.abun), mean(u.abun), sd(b.abun), sd(u.abun), length(b.abun), length(u.abun))
mes(mean(b.rich), mean(u.rich), sd(b.rich), sd(u.rich), length(b.rich), length(u.rich))

# burn + harvest #
mes(35.11, 54.33, 27.84, 29.71, 9, 9) # abun 5-burn
mes(27.89, 54.33, 11.95, 29.71, 9, 9) # abun 30-burn
mes(31.78, 54.33, 14.32, 29.71, 9, 9) # abun 60-burn

mes(6.67, 6.22, 1.66, 2.11, 9, 9) # rich 5-burn
mes(8.44, 6.22, 1.94, 2.11, 9, 9) # rich 30-burn
mes(7.67, 6.22, 1.32, 2.11, 9, 9) # rich 60-burn

##########################
## Colby (Dissertation) ##
##########################
# Carabidae #
# Treatment 1X NA #
mes(12, 2, ) # abundance
mes() # richness

# Ants #
mes(14.19, 4.73, 25.52, 6.43, 19, 9) # Abundance
mes(2.5, 3.5, 1, 2.12, 2, 4) # Richness

#############
## Apigian ##
#############
mes(14.82, -4.65, 2.43, 13.58, 3, 3) # richness
mes(-4.13, -1.43, 3.47, 1.29, 3, 3) # abundance

#################
## Chen et al. ##
#################
# Carabidae #
# Wildfire #
mes(8.19, 1.65, 11.65-8.19, 2.74-1.65, 160, 160) # richness from data thief

# Fire + Mechnanical #
mes(2.32, 1.65, 3.82-2.32, 2.74-1.65, 160, 160) # richness from data thief

#####################
## Roughley et al. ##
#####################

# First Year #
mes(192.25, 139.5, 110.65, 65.39, 4, 4) # spring abundance
mes(282, 139.5, 170.56, 65.39, 4, 4) # summer abundance
mes(137.5, 139.5, 74.04, 65.39, 4, 4) # fall abundance

mes(20.75, 17.75, 1.71, 3.95, 4, 4) # spring richness
mes(20.5, 17.75, 4.36, 3.95, 4, 4) # summer richness
mes(18, 17.75, 5.23, 3.95, 4, 4) # fall richness

# Second Year #
mes(216.5, 171, 110.65, 106.52, 4, 4) # spring abundance
mes(189.75, 171, 170.56, 106.52, 4, 4) # summer abundance
mes(169.75, 171, 74.04, 106.52, 4, 4) # fall abundance

mes(20.5, 21, 5.07, 5.66, 4, 4) # spring richness
mes(24.75, 21, 2.63, 5.66, 4, 4) # summer richness
mes(20.5, 21, 5.92, 5.66, 4, 4) # fall richness

# Third Year #
mes(183.5, 166.5, 46.26, 60.42, 4, 4) # spring abundance
mes(128.25, 166.5, 78.56, 60.42, 4, 4) # summer abundance
mes(103.75, 166.5, 36.2, 60.42, 4, 4) # fall abundance

mes(21.75, 22.5, 2.87, 4.2, 4, 4) # spring richness
mes(20.5, 22.5, 5, 4.2, 4, 4) # summer richness
mes(21.25, 22.5, 3.3, 4.2, 4, 4) # fall richness

# Fourth Year #
mes(127.25, 117.75, 58.67, 39.64, 4, 4) # spring abundance
mes(117.5, 117.75, 56.18, 39.64, 4, 4) # summer abundance
mes(125.75, 117.75, 42.36, 39.64, 4, 4) # fall abundance

mes(18.75, 16.25, 3.77, 5.12, 4, 4) # spring richness
mes(18, 16.25, 5, 3.56, 4, 4) # summer richness
mes(21.25, 16.25, 2.06, 5.12, 4, 4) # fall richness

###########
## Short ##
###########
## Carabidae ##
mes(9.8, 23.6, 14.7/sqrt(3), 32.3/sqrt(3), 3, 3) # abundance 1y
mes(13.5, 15.7, 11.1/sqrt(3), 22.6/sqrt(3), 3, 3) # abundance 2y

## Bees ##
mes(18.4, 16.8, 20.4/sqrt(3), 23.2/sqrt(3), 3, 3) # abundance y1
mes(22.5, 17.3, 17.2/sqrt(3), 17.1/sqrt(3), 3, 3) # abundance y2

## Ants ##
mes(142.4, 92.0, 80.3*sqrt(3), 42.1*sqrt(3), 3, 3) # abundance 1y
mes(81.8, 91.7, 14.8*sqrt(3), 24.7*sqrt(3), 3, 3) # abundance 2y

#############################
## Nelson (Masters Thesis) ##
#############################
# Carabidae #
# Fire #
mes(184.33, 174.56, 57.07, 131.68, 3, 9) # 2017 abundance
mes(16, 15.78, 1.73, 2.82, 3, 9) # 2017 richness

mes(164.33, 93.25, 75.27, 62.39, 6, 4) # 2018 abundance 1y
mes(14.67, 12, 2.16, 6.88, 6, 4) # 2018 richness 1y

mes(146.67, 93.25, 40.38, 62.39, 3, 4) # 2018 abundance 2y
mes(17.67, 12, 1.15, 6.88, 3, 4) # 2018 richness 2y

# Fire + Grazing #
mes(397, 174.56, 392.3869, 131.68, 5, 9) # 2017 abundance y0
mes(17.8, 15.78, 4.55, 2.82, 5, 9) # 2017 richness y0

mes(407.67, 93.25, 96.69, 62.39, 3, 4) # 2018 abundance y0
mes(13.67, 12, 4.16, 6.88, 3, 4) # 2018 richness y0

mes(130.5, 93.25, 20.51, 62.39, 2, 4) # 2018 abundance fire freq y0
mes(14, 12, 1.41, 6.88, 2, 4) # 2018 richness fire freq y0

##########################
## Campbell et al. 2018 ##
##########################
# Carabids #
# Burn + Mechanical #
csd <- sum(c(0.33*sqrt(24), 2.03*sqrt(24))^2)
bsd <- sum(c(1.86*sqrt(24), 1.2*sqrt(24))^2)

mes(3.67+1.67, 0.67+5.67, sqrt(bsd/2), sqrt(csd/2), 24, 24)

# Burn #
c.abun <- c(0.8, 8.0)
b.abun <- c(5.1, 10.0)

c.abun.sd <- c((0.2*sqrt(24))^2, (2.5*sqrt(24))^2)
b.abun.sd <- c((2.0*sqrt(24))^2, (3.3*sqrt(24))^2)
mes(mean(b.abun), mean(c.abun), sum(b.abun.sd)/2, sum(c.abun.sd)/2, 24, 24)

#############################
## Davis (Grey-literature) ##
#############################
mes(18, 18, 10.83, 9.45, 10, 9) # abun low
mes(18, 18, 7.96, 9.45, 10, 9) # abun hi
mes(8, 6.78, 2.54, 2.54, 10, 9) # rich low
mes(9.3, 6.78, 3.30, 2.54, 10, 9) # rich hi

####################
## Neeman & Dafni ##
####################

# n = 120 days #
mes(mean(1.09,0), mean(0.47, 0.63), sqrt((2.7^2+.007^2)/(120)), sqrt((0.95^2+1.1^2)/(120)), 120, 120) # abundance

####################
## Dodge (Thesis) ##
####################

abun.b <- mean(c(6.5, 5.33, 2.86, 12.67, 20.33, 17))
abun.ub <- mean(c(3.14, 3.14, 0.29, 17.57, 13.71, 3))
rich.b <- mean(c(4.83, 4.5, 2.43, 6.17, 7.33, 4.57))
rich.ub <- mean(c(2.86, 3, 0.57, 5.29, 5.29, 1.57))

abun.se.b <- sum(c(1.23, 1.86, 1.01, 5.22, 4.2, 3.47))/7
abun.se.ub <- sum(c(1.53, 0.99, 0.18, 6.51, 4.92, 1.05))/7
rich.se.b <- sum(c(0.31, 0.67, 0.75, 1.64, 0.71, 0.75))/7
rich.se.ub <- sum(c(1.03, 0.69, 0.37, 1.61, 1.25, 0.43))/7

mes(abun.b, abun.ub, abun.se.b, abun.se.ub, 6, 6)
mes(rich.b, rich.ub, rich.se.b, rich.se.ub, 6, 6)

#######################
## Welti et al. 2018 ##
#######################
## Bees ##
## Fire ##
mes(22, 8.5, 5.7, 2.1, 4, 2) #1y rich
mes(269, 54, 11, 26, 4, 2) #1y abun
mes(25.5, 8.5, 3.5, 2.1, 4, 2) #4y rich
mes(187, 54, 53, 26, 4, 2) #4y abun

## Fire + Grazing ##
mes(10.25, 9, 4.03, 1.41, 4, 2) # 1y rich
mes(78.5, 53.5, 20.42, 24.75, 4, 2) # 1y abun
mes(11.5, 9, 3.697, 1.41, 4, 2) # 4y rich
mes(88.5, 53.5, 43.26, 24.75, 4, 2) # 4y abun

## Ants ##
## Fire + Grazing ##
mes(1.25, 2, 1.5, 0, 4, 2) # 1y rich
mes(3.25, 13.5, 3.95, 7.78, 4, 2) # 1y abun
mes(2.25, 3, 0.95, 1.41, 4, 2) # 4y rich
mes(8, 9.5, 6.38, 10.61, 4, 2) # 4y abun

## Butterflies ##
## Fire ##
mes(4, 3, 0, 1.4, 4, 2) #1y rich
mes(10, 9.5, 8.5, 11, 4, 2) #1y abun
mes(3.5, 3, 3.5, 1.4, 4, 2) #4y rich
mes(6.5, 9.5, 7.8, 11, 4, 2) #4y abun

## Fire + Grazing ##
mes(0.75, 3, 0.5, 1.41, 4, 2) # 1y rich
mes(1.5, 9.5, 1.73, 10.61, 4, 2) # 1y abun
mes(2.25, 0.96, 3, 1.41, 4, 2) # 4y rich
mes(8, 6.38, 9.5, 10.61, 4, 2) # 4y abun

#########################
## Carbone et al. 2017 ##
#########################

## ADD FORMULA FOR SE TO SD!!!!!! ##
mes(4.9, 4.4, (5.5-4.3)/3.92, (6.2-2.6)/3.92, ) #ub vs low y1 s1
mes(2.7, 4.4, (3.2-2.3)/3.92, (6.2-2.6)/3.92, ) #ub vs hi y1 s1

mes(5.9, 3.8, (9.8-2.1)/3.92, (5.1-2.3)/3.92, ) #ub vs hi y2 s1
mes(1.5, 3.8, (1.8-1.3)/3.92, (5.1-2.3)/3.92, ) #ub vs hi y2 s1

mes(3.5, 3.9, (4-3)/3.92, (4.6-3.2)/3.92, ) #ub vs low y1 s2
mes(3.5, 3.9, (4-3)/3.92, (4.6-3.2)/3.92, ) #ub vs hi y1 s2

mes(1.7, 3.7, (1.8-1.5)/3.92, (5.3-2)/3.92, ) #ub vs hi y2 s2
mes(5, 3.7, (5.9-3.2)/3.92, (5.3-2)/3.92, ) #ub vs hi y2 s2

###################
## Nuland et al. ##
###################

mes(3.6, 2.1, 2.7, 2.12, 30, 30) # abundance
mes(1.1, 0.97, 0.66, 0.76, 30, 30) # richness

#########################
## Lettow et al. 2018  ##
#########################
## Burn Data in Table 1 ##
mes(41, 37, 5, 3, 20, 20) # 2011 abundance
mes(11, 12, 2, 1, 20, 20) # 2011 richness
mes(38, 31, 1, 5, 20, 20) # 2012 abundance
mes(13, 12, 2, 1, 20, 20) # 2012 richness

## Mechanical Burn Data in Table 1 ##
mes(48, 37, 9, 3, 20, 20) # 2011 abundance
mes(17, 12, 2, 1, 20, 20) # 2011 richness
mes(54, 31, 8, 5, 20, 20) # 2012 abundance
mes(21, 12, 3, 1, 20, 20) # 2012 richness

######################
## Love et al. 2016 ##
######################
burn <- c(15, 58, 12, 21) # abundance
control <- c(52, 47, 24, 27) # abundance
mes(mean(burn), mean(control), sd(burn), sd(control), 4, 4)

burn <- c(4, 8, 5, 5) # richness
control <- c(10, 19, 9, 4) # richness
mes(mean(burn), mean(control), sd(burn), sd(control), 4, 4)

##########################
## LoPresti et al. 2018 ##
##########################
controlsN <-c(5.76, 3.86, 3.33, 2.60, 2.40, 2.25, 2.11, 2.00, 2.00, 1.47)
controlsSD <-c(5.83, 6.55, 2.08, 1.14, 1.34, 1.91, 1.47, 0.0, 1.41, 1.06)

PreburnN <-c(5.36, 2.00, 1.18, 1.61, 1.89, 2.00, 1.58, 1.52, 1.00, 1.00, 1.00)
PreburnSD <-c(4.62, 0.0, 0.53, 1.02, 1.41, 1.67, 0.99, 0.85, 0.0, 0.0, 0.0)

controlsM<-cbind(controlsN, controlsSD)
PreburnM<-cbind(PreburnN, PreburnSD)

controlsPSD<-sqrt(
  ((controlsM[1,1] - 1)*controlsM[1,2]^2 + 
     (controlsM[2,1] - 1)*controlsM[2,2]^2 +
     (controlsM[3,1] - 1)*controlsM[3,2]^2 +
     (controlsM[4,1] - 1)*controlsM[4,2]^2 +
     (controlsM[5,1] - 1)*controlsM[5,2]^2 +
     (controlsM[6,1] - 1)*controlsM[6,2]^2 +
     (controlsM[7,1] - 1)*controlsM[7,2]^2 +
     (controlsM[8,1] - 1)*controlsM[8,2]^2 +
     (controlsM[9,1] - 1)*controlsM[9,2]^2 +
     (controlsM[10,1] - 1)*controlsM[10,2]^2)/
    (sum(controlsM[,1])-nrow(controlsM))
)
                  
PreburnPSD<-sqrt(
  abs(((PreburnM[1,1] - 1)*PreburnM[1,2]^2 + 
         (PreburnM[2,1] - 1)*PreburnM[2,2]^2 +
         (PreburnM[3,1] - 1)*PreburnM[3,2]^2 +
         (PreburnM[4,1] - 1)*PreburnM[4,2]^2 +
         (PreburnM[5,1] - 1)*PreburnM[5,2]^2 +
         (PreburnM[6,1] - 1)*PreburnM[6,2]^2 +
         (PreburnM[7,1] - 1)*PreburnM[7,2]^2 +
         (PreburnM[8,1] - 1)*PreburnM[8,2]^2 +
         (PreburnM[9,1] - 1)*PreburnM[9,2]^2 +
         (PreburnM[10,1] - 1)*PreburnM[10,2]^2 +
         (PreburnM[11,1] - 1)*PreburnM[11,2]^2)/
        (sum(PreburnM[,1])-nrow(PreburnM)))
)

mes(sum(controlsM[,1]), sum(PreburnM[,1]), controlsPSD*sqrt(10), PreburnPSD*sqrt(11), 10, 11)

#############################
## Potts et al. 2006       ##
#############################
mes(18, 17, 2.4*sqrt(3), 2.0*sqrt(3), 3, 3) # richness
mes(469, 161, 70*sqrt(3), 32*sqrt(3), 3, 3) # abundance

mes(15, 17, 2.3*sqrt(3), 2.0*sqrt(3), 3, 3) # richness
mes(354, 161, 81*sqrt(3), 32*sqrt(3), 3, 3) # abundance

#############################
## Peralta et al.          ##
#############################
fb <- c(6, 4) # richness
ib <- c(3, 4, 9, 7) # richness
ob <- c(4, 3, 5, 7) # richness
control <- c(7, 5, 5, 7) # richness

mes(mean(fb), mean(control), sd(fb), sd(control), length(fb), length(control))
mes(mean(ib), mean(control), sd(ib), sd(control), length(ib), length(control))
mes(mean(ob), mean(control), sd(ob), sd(control), length(ob), length(control))

#############################
## Moretti et al.          ##
#############################
mes(20, 12.2, 1.9*sqrt(6), 1.1*sqrt(6), 6, 6) # richness med
mes(64.8, 45.3, 7.1*sqrt(6), 6.1*sqrt(6), 6, 6) # abundance med

mes(11.5, 12.2, 2.3*sqrt(6), 0.7*sqrt(6), 6, 6) # richness med in
mes(40, 45.3, 6.3*sqrt(6), 6.1*sqrt(6), 6, 6) # abundance med in

mes(9.7, 12.2, 2*sqrt(3), 1.1*sqrt(6), 3, 6) # richness med ma
mes(39, 45.3, 40*sqrt(3), 6.1*sqrt(6), 6, 6) # abundance med ma


mes(67, 48.8, 6.4*sqrt(5), 1.3*sqrt(6), 5, 6) # richness temp
mes(928, 504, 234*sqrt(5), 63*sqrt(6), 5, 6) # abundance temp

mes(52,  48.8, 1.8*sqrt(7), 1.3*sqrt(6), 7, 6) # richness temp in
mes(554, 504, 103*sqrt(7), 63*sqrt(6), 7, 6) # abundance temp in

mes(49.3,  48.8, 3.4*sqrt(3), 1.3*sqrt(6), 3, 6) # richness temp ma
mes(508, 504, 158*sqrt(3), 63*sqrt(6), 3, 6) # abundance temp ma


#############################
## DISSERTATION: Simanonok ##
#############################
# Bees #
mes(1.83, 0, 0.75, 0, 6, 6) # 3y Hi Rich
mes(1.61, 0, 1.30, 0, 6, 6) # 3y Hi Abun

mes(1.83, 0, 1.17, 0, 6, 6) # 3y Mx Rich
mes(3.09, 0, 1.35, 0, 6, 6) # 3y Mx Abun

mes(2.5, 0, 1.38, 0, 6, 6) # 4y Hi Rich
mes(1.1, 0, 0.37, 0, 6, 6) # 4y Hi Abun

mes(1.5, 0, 1.7, 0, 6, 6) # 4y Mx Rich
mes(0.86, 0, 0.87, 0, 6, 6) # 4y Mx Abun

mes(1.8, 0, 1.72, 0, 6, 6) # 25y Hi Rich
mes(0.696, 0, 0.75, 0, 6, 6) # 25y Hi Abun

mes(0.6, 0, 0.54, 0, 6, 6) # 25y Mx Rich
mes(0.14, 0, 0.17, 0, 6, 6) # 25y Mx Abun

mes(3, 0, 0.63, 0, 6, 6) # 9y Hi Rich
mes(0.78, 0, 0.37, 0, 6, 6) # 9y Hi Abun

mes(3, 0, 1.09, 0, 6, 6) # 9y Mx Rich
mes(1.47, 0, 0.41, 0, 6, 6) # 9y Mx Abun

####################
## Granath et al. ##
####################
## Data on Dryad ##
control <- c(0,0,0,0,2,2,6,0,9,1,2,2) 
burn <- c(7,7,3,2,2,11,6,2,5,6,4,3)
mes(mean(burn), mean(control), sd(burn), sd(control), 12, 12) # abundance

control <- c(0,0,0,0,2,2,2,0,4,1,3,2)
burn <- c(7,6,3,1,2,6,4,2,4,3,2,2)
mes(mean(burn), mean(control), sd(burn), sd(control), 12, 12) # richness

## Fire + Mechnical Burn ##
mes(12, 3.17, 3.28, 1.64, 12, 12) # 0-ret richness
mes(35, 5.58, 15.93, 3.579, 12, 12) # 0-ret abundance

mes(13, 3.17, 2.73, 1.64, 12, 12) # 10-ret richness
mes(34.92, 5.58, 13.86, 3.579, 12, 12) # 10-ret abundance

mes(13.33, 3.17, 3.92, 1.64, 12, 12) # 50-ret richness
mes(61.75, 5.58, 59.29, 3.579, 12, 12) # 50-ret abundance

###################
## Burkle et al. ##
###################
# inverse effect

b.abun <- c(2.7, 1.2, 1.8)
b.abunSD <- c(0.2, 0.1, 0.1)
u.abun <- c(5.4, 3.1, 1.9)
u.abunSD <- c(0.2, 0.2, 0.1)

mes(mean(b.abun), mean(u.abun), sqrt((sum(b.abunSD))/3), sqrt((sum(u.abunSD))/3), 36, 17)

b.rich <- c(2.1, 1.2, 1.6)
b.richSD <- c(0.2, 0.03, 0.1)
u.rich <- c(4.0, 2.6, 1.8)
u.richSD <- c(0.1, 0.1, 0.03)

mes(mean(b.rich), mean(u.rich), sqrt((sum(b.richSD))/3), sqrt((sum(u.richSD))/3), 36, 17)

###################
## Nuland et al. ##
###################

# this paper is fire freq!

mes(1.133, 1, 0.681, 0.871, 30, 30) # richness
mes(3.633, 2.167, 2.735, 2.198, 30, 30) # abundance

###########
## Pryke ##
###########

cbb <- c(1,16,0)
cbu <- c(0,18,135,0)
fb <- c(1,14,30)
fu <- c(0,36,36,49)

mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 4)
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

cbb <- c(1,8,0)
cbu <- c(0,6,12,0)
fb <- c(1,6,7)
fu <- c(0,8,8,3)

mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 4)
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

######################
## Tucker and Rehan ##
######################

mes(22.69, ) # abundance
mes() # richness

####################
## Adedoja et al. ##
####################
mes(14, 24.33, 13, 18.90326, 3, 3) # abundance hilltop
mes(12, 31.667, 8.544, 27.68273, 3, 3) # abundance hill slope
mes(15, 22, 10.54, 12.28821, 3, 3) # abundance valley

mes(3, 4.667, 2.65, 4.04, 3, 3) # richness hilltop
mes(2.667, 3, 1.15, 1, 3, 3) # richness hill slope
mes(3, 2.667, 1, 1.53, 3, 3) # richness valley

####################
## Bogusch et al. ##
####################
mes(8.87, 4.75, 7.9, 3.86, 4, 15) # abundance
mes(5, 2.75, 3.12, 1.5, 4, 15) # richness

## BUTTERFLIES ##

##################
## Carbone 2017 ##
##################
mes() # richness 

#######################
## Welti et al. 2018 ##
#######################


###################
## Nuland et al. ##
###################
mes(0.03, 0, 0.18, 0, 30, 30)

#########################
## Campbell et al 2007 ##
#########################
## Bees ##
bAbun <- c(104.3, 4.7, 5.0)
bAbunVar <- c(25.2, 2.7, 2.0)
bAbunSD <- bAbunVar * sqrt(10)
uAbun <- c(118, 4.0, 5.0)
uAbunVar <- c(52.9, 2.0, 1.0)
uAbunSD <- uAbunVar * sqrt(10)
mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10) # year 1

bAbun <- c(152.3, 5.7, 2.0)
bAbunVar <- c(73, 0.88, 1.0)
bAbunSD <- bAbunVar * sqrt(10)
uAbun <- c(106.7, 17.7, 3.0)
uAbunVar <- c(44.2, 11.2, 1.0)
uAbunSD <- uAbunVar * sqrt(10)
mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10) # year 2

## Mechanical Burn Table 2 ##
bAbun <- c(238.7, 14.3, 16)
bAbunVar <- c(37, 4.5, 1.2)
bAbunSD <- bAbunVar * sqrt(10)
uAbun <- c(118, 4.0, 5.0)
uAbunVar <- c(52.9, 2.0, 1.0)
uAbunSD <- uAbunVar * sqrt(10)
mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10) # year 1

bAbun <- c(341.7, 14.7, 13.3)
bAbunVar <- c(96.4, 1.2, 4.5)
bAbunSD <- bAbunVar * sqrt(10)
uAbun <- c(106.7, 17.7, 3.0)
uAbunVar <- c(44.2, 11.2, 1.0)
uAbunSD <- uAbunVar * sqrt(10)
mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10) # year 2

## Butterflies ##
## Data in table 2 ##
bAbun <- c(13.7, 3.5)
bAbunVar <- c(3.5, 0.5)
bAbunSD <- bAbunVar * sqrt(10)
uAbun <- c(11.3, 5.3)
uAbunVar <- c(4.9, 1.5)
uAbunSD <- uAbunVar * sqrt(10)
mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10) # year 1

bAbun <- c(20.7, 4.0)
bAbunVar <- c(10.7, 2.5)
bAbunSD <- bAbunVar * sqrt(10)
uAbun <- c(11.3, 4.0)
uAbunVar <- c(2.7, 0)
uAbunSD <- uAbunVar * sqrt(10)
mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10) # year 2

## Mechanical Burn Table 2 ##
bAbun <- c(25, 8.5)
bAbunVar <- c(6.1, 6.5)
bAbunSD <- bAbunVar * sqrt(10)
uAbun <- c(11.3, 5.3)
uAbunVar <- c(4.9, 1.5)
uAbunSD <- uAbunVar * sqrt(10)
mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10) # year 1

bAbun <- c(20.7, 10.3)
bAbunVar <- c(6.8, 3.0)
bAbunSD <- bAbunVar * sqrt(10)
uAbun <- c(11.3, 4.0)
uAbunVar <- c(2.7, 0)
uAbunSD <- uAbunVar * sqrt(10)
mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10) # year 2

#####################
## Fleishman et al ##
#####################
# butterflies #
bRich <- c(18, 28, 28, 20, 15, 22, 16)
uRich <- c(43, 47, 54, 38, 46)

mes(mean(bRich), mean(uRich), sd(bRich), sd(uRich), 7, 5)

###################
## Warchola 2015 ##
###################

# abundance
mes(0.048, 0.030, sqrt(80)*(0.068-0.034)/3.92, sqrt(80)*(0.055-0.016)/3.92, 80, 80)

###################
## Warchola 2015 ##
###################

# Abundance
mes(9.74, 7.48, 7.88, 6.69, 67, 91)
mes(24.99, 18.72, 22.41, 22.89, 67, 91)
mes(15.57, 16.36, 18.44, 19.74, 67, 91)

###########
## Pryke ##
###########

cbb <- c(1,4,7)
cbu <- c(0,4,8,3)
fb <- c(0,6,1)
fu <- c(0,0,0,4)

mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 4)
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

cbb <- c(3,5,3)
cbu <- c(4,5,5,5)
fb <- c(5,5,3)
fu <- c(4,5,5,5)

mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 4)
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

#################
## Kwon et al. ##
#################

# Richness #
mes(27, 30, sqrt(6)*4, sqrt(6)*3, 6, 6) # 0 year
mes(63, 45, sqrt(7)*7, sqrt(7)*6, 7, 7) # 1 year
mes(30, 46, sqrt(7)*5, sqrt(7)*6, 7, 7) # 2 year
mes(37, 46, sqrt(6)*4, sqrt(6)*2, 6, 6) # 3 year
mes(41, 41, sqrt(5)*5, sqrt(5)*6, 5, 5) # 4 year

###########################
## Mciver and Macke 2014 ##
###########################

## year one
mes(18.4, 5.4, 29.72, 9.26, 6, 6) # JP Abun
mes(10.32, 14.4, 12.17, 5.67, 5, 6) # PJ Abun
mes(9.2, 3.6, 5.67, 1.2, 3, 3) # SBE Abun
mes(13, 11.8, 9.94, 7.14, 6, 6) # SBW Abun
mes(13.44, 12.48, 19.63, 13.13, 5, 5) # WJ Abun

mes(2.67, 1.83, 2.34, 2.23, 6, 6) # JP Rich
mes(2, 4, 1.41, 2.97, 5, 6) # PJ Rich
mes(3.33, 1.33, 3.21, 0.58, 3, 3) # SBE Rich
mes(2, 1.83, 0.89, 0.41, 6, 6) # SBW Rich
mes(3, 3, 1.58, 1.41, 5, 5) # WJ Rich

## year two
mes(18.24, 8.16, 10.44, 6.08, 5, 5) # JP Abun
mes(11.53, 28.08, 11.93, 23.58, 6, 5) # PJ Abun
mes(8.1, 8.4, 5.49, 4.596, 4, 4) # SBE Abun
mes(10.15, 14.8, 4.3, 18.24, 6, 6) # SBW Abun
mes(27.78, 22.32, 32.2, 31.53, 5, 5) # WJ Abun

mes(2.8, 1.6, 1.48, 1.52, 5, 5) # JP Rich
mes(3.5, 2.8, 2.43, 1.64, 6, 5) # PJ Rich
mes(2.25, 2.5, 0.96, 1.91, 4, 4) # SBE Rich
mes(3.17, 3.17, 1.72, 2.64, 6, 6) # SBW Rich
mes(4.8, 3.2, 2.77, 2.49, 5, 5) # WJ Rich

## year three
mes(18.69, 8.57, 15.53, 13.92, 7, 7) # JP Abun
mes(24.6, 36.8, 21.98, 33.47, 6, 6) # PJ Abun
mes(7.15, 3, 2.94, 5.23, 4, 4) # SBE Abun
mes(6.25, 2.8, 3.50, 2.36, 7, 6) # SBW Abun
mes(32.23, 75.43, 45.23, 107.5, 7, 7) # WJ Abun

mes(2.86, 1.57, 2.67, 2.299, 7, 7) # JP Rich
mes(4.17, 4.33, 1.72, 2.25, 6, 6) # PJ Rich
mes(3.25, 1.25, 1.5, 1.89, 4, 4) # SBE Rich
mes(1.86, 1.5, 0.69, 1.38, 7, 6) # SBW Rich
mes(4, 4.14, 2.65, 3.13, 7, 7) # WJ Rich

## year four
mes(7.2, 5.2, 2.49, 4.34, 6, 6) # JP Abun
mes(24.96, 20.88, 13.82, 18.41, 5, 5) # PJ Abun
mes(5.28, 1.44, 4.21, 1.13, 5, 5) # SBE Abun
mes(9.6, 5.14, 8.28, 4.79, 6, 7) # SBW Abun
mes(16.29, 65.1, 22.03, 127.15, 7, 8) # WJ Abun

mes(3.17, 2.33, 1.33, 2.25, 6, 6) # JP Rich
mes(5.4, 3.2, 2.88, 1.79, 5, 5) # PJ Rich
mes(2.4, 0.6, 1.5, 0.55, 5, 5) # SBE Rich
mes(2.67, 2.14, 2.88, 1.68, 6, 7) # SBW Rich
mes(3.29, 3.63, 1.98, 2.72, 7, 8) # WJ Rich

## year five
mes(21.84, 5.04, 32.94, 4.75, 5, 5) # JP Abun
mes(17.1, 24, 12.98, 20.39, 4, 4) # PJ Abun
mes(1.6, 2.8, 2.77, 3.02, 3, 3) # SBE Abun
mes(4.8, 6.72, 4.07, 8.17, 5, 5) # SBW Abun
mes(16.8, 22.4, 13.34, 16.48, 6, 6) # WJ Abun

mes(3.8, 2.6, 4.66, 1.82, 5, 5) # JP Rich
mes(3.25, 3, 2.87, 1.83, 4, 4) # PJ Rich
mes(1, 1.33, 1.73, 1.15, 3, 3) # SBE Rich
mes(1.6, 1.8, 1.52, 1.92, 5, 5) # SBW Rich
mes(3.33, 3.5, 1.97, 1.22, 6, 6) # WJ Rich

## year six
mes(16.5, 15.9, 10.79, 15.03, 4, 4) # JP Abun
mes(11.1, 14.7, 3.71, 11.06, 4, 4) # PJ Abun
mes(6, 3.2, 1.2, 3.67, 3, 3) # SBE Abun
mes(30.6, 20.1, 17.35, 3.96) # WJ Abun

mes(3.5, 2.75, 1.73, 1.71, 4, 4) # JP Rich
mes(2.5, 2.75, 1.29, 2.36, 4, 4) # PJ Rich
mes(2, 1, 1.73, 1, 3, 3) # SBE Rich
mes(4.5, 4, 1.29, 1.41, 4, 4) # WJ Rich

## year seven
mes(4.8, 4.8, 4.49, 1.96, 4, 4) # JP Abun
mes(51.6, 50.4, 6.79, 11.88, 2, 2) # PJ Abun
mes(3, 1.2, 0.85, 1.697, 2, 2) # SBE Abun
mes(15.12, 14.64, 15.57, 17.17, 5, 5) # WJ Abun

mes(2.5, 2.75, 2.08, 1.26, 4, 4) # JP Rich
mes(5.5, 2.5, 2.12, 0.71, 2, 2) # PJ Rich
mes(3.2, 1.4, 2.49, 1.14, 5, 5) # WJ Rich

######################
## Riva et al, 2018 ##
######################

mes(0.005556, 0.358333, 0.074536, 0.683489, 180, 120) 

######################
## Scandurra et al. ##
######################

mes(15.3, 12.8, 0.08, 1, 24, 24) # richness
mes(131.5, 71.1, 14.6, 9.8, 24, 24) # abundance

###################
## Serrat et al. ##
###################

mes(56.9, 74.3, 42.9, 20.5, 70, 70) # Abudance
mes(16.0, 13.9, 5, 8.2, 70, 70) # Richness

###############
## Fleishman ##
###############

mes(16.6, 32.6, 6.025, 6.5803, 5, 5) # 2007
mes(17, 35.8, 4.42, 7.5299, 5, 5) # 2008

#################################
## Cleary & Mooers 2004 & 2006 ##
#################################

mes(28.76, 53.91, sqrt(3200)*(8.22)/3.92, sqrt(2841)*(5.35)/3.92, 3200, 2841) # species abundance from 2006
esc_t(t=-1.567, grp1n=3, grp2n=3, es.type="g") # richness from 2006

mes(64.47, 93.49, sqrt(3200)*(103.86-85.3)/3.92, sqrt(2841)*(75.67-60.48)/3.92, 3200, 2841) # species richness from 2004

############
## Powell ##
############
mes(0.00872, .0322, 0.035, 0.035, 62, 19) # abundance


#############
## Nowicki ##
#############
mes(2.20, 1.84, 1.98, 1.48, 14, 41) # abundance

############
## Dollar ##
############
mes(5.69, 6.41, sqrt(6)*(6.94-5.69), sqrt(6)*(7.66-6.41), 6, 6) # abundance 2008
mes(17.0, 19.74, sqrt(6)*(18.97-17.0), sqrt(6)*(21.12-19.74), 6, 6) # richness 2008
mes(9.32, 7.93, sqrt(6)*(10.87-9.32), sqrt(6)*(9.45-7.93), 6, 6) # abundance 2009
mes(14.85, 12.29, sqrt(6)*(16.72-14.85), sqrt(6)*(14.17-12.29), 6, 6) # richness 2009

#################
# Bohls et al. ##
#################
# Butterflies #
mes(5.27, 20.83, 4*1.35, 4*1.24, 16, 16) # abundance
mes(3.00, 5.00, 4*0.58, 4*0.00, 16, 16) # richness

###################
## Beaumont 2012 ##
###################
# ants #
# data in table 1
mes(651.3, 291.2, 89.81*5, 48.13*5, 25, 25) # 2006 abundance
mes(537.0, 165.5, 107.87*5, 28.22*5, 25, 25) # 2007 abundance
mes(418.3, 187.5, 55.02*5, 37.4*5, 25, 25) # 2008 abundance

mes(11.9, 13.1, 1.25*5, 0.9*5, 25, 25) # 2006 richness
mes(11.1, 13.0, 0.57*5, 0.96*5, 25, 25) # 2007 richness
mes(10.8, 12.6, 1.36*5, 1.36*5, 25, 25) # 2008 richness

##############
## Atchison ##
##############
esc_f(f=14.4, grp1n=30, grp2n=30, es.type="g")

###########
## Anjos ##
###########
esc_f(f=11.1, grp1n=144, grp2n=144, es.type="g") # abundance
esc_f(f=0.6, grp1n=144, grp2n=144, es.type="g") # richness

############
## Andrew ##
############
esc_f(f=0.39, grp1n=36, grp2n=36, es.type="g") # richness
esc_f(f=4.61, grp1n=36, grp2n=36, es.type="g") # abundance

#############
## Apigian ##
#############
mes(531.33, 756.44, 572.62, 1189, 20, 20) # abundance

###################
## Sackmann 2006 ##
###################
c.abun <- c(16.4, 1, 6.5, 9.3, 9, 0.5, 36.4, 0, 1, 0.87, 3, 6, 3, 1)
b.abun <- c(41.2, 42.27, 4.3, 1.3, 41.8, 0, 0.5, 205, 16.6, 0, 15.6, 14.5, 0.6, 4, 0)
c.abun.sd <- c(23.2, 1, 10.6, 14.5, 6.5, 0, 0.5, 68.2, 0, 0, 1.3, 4.2, 1.7, 4.2, 0)^2
b.abun.sd <- c(25.4, 0, 7.2, 2.3, 42, 0, 0.5, 166.1, 26.2, 0, 19, 19.1, 0.8, 4.2, 0)^2
mes(mean(b.abun), mean(c.abun), sqrt(sum(b.abun.sd)/8), sqrt(sum(c.abun.sd)/8), 8, 8) # abundance forest

c.abun <- c(28.4, 1, 17.5, 4, 0, 2.0, 539, 0, 2, 1, 43, 0, 1)
b.abun <- c(24.8, 0, 11.3, 0.3, 1, 1, 773, 6722, 656, 25.8, 0, 1, 0)
c.abun.sd <- c(31.4, 0, 13.4, 3.6, 0, 1.4, 1079, 0, 1.8, 1.7, 0, 0, 0)
b.abun.sd <- c(20, 0, 5, 0.6, 0, 0, 1433, 0, 1250, 25, 0, 0, 0)
mes(mean(b.abun), mean(c.abun), sqrt(sum(b.abun.sd)/6), sqrt(sum(c.abun.sd)/6), 6, 6) # abundance scrub

c.abun <- c(25.7, 4.0, 43, 11318, 0, 458, 10.5, 11.4, 0, 11)
b.abun <- c(38.5, 2.7, 6, 1, 7238, 79.3, 365, 28, 14.8, 26, 0)
c.abun.sd <- c(49.5, 6.9, 0, 43, 22265, 0, 785, 14.8, 0, 11.2, 0, 0)
b.abun.sd <- c(75.7, 1.5, 0, 0, 11640, 121.8, 601, 28.3, 0, 20.6, 0, 0)
mes(mean(b.abun), mean(c.abun), sqrt(sum(b.abun.sd)/5), sqrt(sum(c.abun.sd)/5), 5, 5) # abundance steppe

###################
## Hanula et al. ##
###################
c.abun <- c(726.5, 803.17)
b.abun <- c(978.7, 764.4)

c.abun.sd <- c((69.04*sqrt(6))^2, (54.81*sqrt(6))^2)
b.abun.sd <- c((75.3*sqrt(6))^2, (37.5*sqrt(6))^2)

mes(mean(b.abun), mean(c.abun), sum(b.abun.sd)/12, sum(c.abun.sd)/12, 12, 12)

###########
## Pryke ##
###########
# y1
cbb <- c(63, 1.69, 1.27, 0)
cbu <- c(41, 1.2, 1.77, 0)
fb <- c()
fu <- c()

mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 4, 4) # abun
mes(mean(fb), mean(fu), sd(fb), sd(fu), 4, 4)

cbb <- c(4.66, 0.4, 0.67)
cbu <- c(2.15, 0.74, 0.92)
fb <- c()
fu <- c()

mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 4, 4) # rich
mes(mean(fb), mean(fu), sd(fb), sd(fu), 4, 4)

#y2
cbb <- c(64.81, 1.85, 0.12)
cbu <- c(51.57, 6.04, 0.33)
mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 3) # abun

cbb <- c(4.95, 0.52, 0.12)
cbu <- c(4.39, 0, 0.04)
mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 3) # rich

#################
## Rice et al. ##
#################
b.abun <- c(0, 75, 40, 20)
c.abun <- c(0, 0, 0, 0)

mes(mean(b.abun), mean(c.abun), sd(b.abun), sd(c.abun), 4, 4) # abundance

#################
## Bess et. al ##
#################

mes(11.95238, 9.208333, 6.01962, 5.999855, 24, 21) # richness
mes(212.0952, 97.41667, 216.1018, 121.361, 24, 21) # abundance

#######################
## Welti et al. 2018 ##
#######################
mes(1.8, 2, 1, 0, 4, 2) #y1 rich
mes(2.5, 13.5, 1.3, 7.78, 4, 2) #y1 abun
mes(2, 2, 0.8, 0, 4, 2) #y4 rich
mes(8, 13.5, 8.4, 7.78, 4, 2) #y4 abun



############
## Barrow ##
############

mes(17.17, 14.33, sqrt(360)*1.85, sqrt(360)*0.80, 360, 360) # dry spin richness
mes(15.83, 19.83, sqrt(360)*0.70, sqrt(360)*1.19, 360, 360) # wet spin richness

mes(22.83, 24.17, sqrt(360)*1.53, sqrt(360)*1.85, 360, 360) # dry sand richness
mes(31.17, 29.58, sqrt(360)*1.66, sqrt(360)*1.01, 360, 360) # wet sand richness

mes(1358.33, 490, sqrt(360)*420.94, sqrt(360)*162.58, 360, 360) # dry spin abundance
mes(3280.83, 894.83, sqrt(360)*324.19, sqrt(360)*235.76, 360, 360) # wet spin abundance

mes(736, 1166.33, sqrt(360)*189.28, sqrt(360)*310.14, 360, 360) # dry sand abundance
mes(1927.83, 2066.06, sqrt(360)*465.41, sqrt(360)*290.34, 360, 360) # wet sand abundance

######################
## Alves-Silva 2011 ##
######################
mes(8.6, 1.2, 4.6, 0.7, 4, 4) # abundance

##############
## Day 2018 ##
##############

mes(8.6, 4.75, (10.4-8.6)*sqrt(5), (5.36-4.75)*sqrt(5), 5, 5) # 2014 abun
mes(15.2, 13.2, (18.4-15.2)*sqrt(5), (15.6-13.2)*sqrt(5), 5, 5) # 2015 abun
mes(22.4, 22.3, (26-22.4)*sqrt(5), (26-22.3)*sqrt(5), 5, 5) # 2016 abun

mes(8.1, 16.8, (10.6-8.1)*sqrt(5), (19.6-16.8)*sqrt(5), 5, 5) # 2014 rich
mes(16.1, 28.6, (18.2-16.1)*sqrt(5), (33.7-28.6)*sqrt(5), 5, 5) # 2015 rich
mes(21.3, 32, (23.9-21.3)*sqrt(5), (36.2-32)*sqrt(5), 5, 5) # 2016 rich

##############
# Orthoptera #
##############

###########
## Pryke ##
###########
#y1
cbb <- c(1.2, 0.2, 0.058)
cbu <- c(1.425, 0.233, 0)
fb <- c()
fu <- c()
mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 3) # abun
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

cbb <- c(1.022, 0.333, 0.058)
cbu <- c(1.15, 0.467, 0)
fb <- c()
fu <- c()
mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 3, 3) # rich
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

#y3
cbb <- c(1.23, 0)
cbu <- c(1.45, 0)
fb <- c()
fu <- c()
mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 2, 2) # abun
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

cbb <- c(1.02, 0)
cbu <- c(1.27, 0)
fb <- c()
fu <- c()
mes(mean(cbb), mean(cbu), sd(cbb), sd(cbu), 2, 2) # rich
mes(mean(fb), mean(fu), sd(fb), sd(fu), 3, 4)

###################
## Hanula et al. ##
###################
c.abun <- c(61.84, 232.33)
b.abun <- c(604.1, 257.9)

c.abun.sd <- c((9.09*sqrt(6))^2, (45.11*sqrt(6))^2)
b.abun.sd <- c((55.82*sqrt(6))^2, (30.54*sqrt(6))^2)

mes(mean(b.abun), mean(c.abun), sum(b.abun.sd)/12, sum(c.abun.sd)/12, 12, 12)

## Bess et. al ##
mes(81.47619, 216.5417, 100.9836, 223.5945, 24, 21) # abudance
mes(4.666667, 5.0, 1.391585, 1.224745, 24, 21) # richness

## Parmenter et. al ##
mes(1.753247, 3.956879, 1.39252, 7.702648, 310, 505) # abundance
mes(0.987138, 1.010121, 0.159869, 0.195885, 310, 505) # richness
