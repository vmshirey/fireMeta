#####################################################################################
## Script for calculating effect sizes from different statistical tests.           ##
## Author: Vaughn M. Shirey                                                        ##
#####################################################################################

require(compute.es)

#################
## Bess et. al ##
#################

#######################################################################################
## Data was lifted from figure 3, on page 780. Control is all pre-fire points; burn  ##
## is all post-fire points (including fire).                                         ##
## Used DataThief, rounded to nearest integer value for richness and hundreths-place ##
## for abundance.                                                                    ##
#######################################################################################

control.abd <- c(32.31, 126.15, 283.08, 282.05, 278.97, 34.87, 9.23, 43.08, 65.64)
control.rch <- c(13, 18, 28, 34, 30, 16, 11, 22, 35)
  
wf.abd <-c(63.08, 102.56, 84.62, 88.21, 50.77, 15.38, 6.67)
wf.rch <- c(29, 31, 22, 20, 18, 12, 20)

print(paste("Control Rch Mean: ", mean(control.rch), " | Control Rch SE: ", sd(control.rch)/sqrt(length(control.rch))))
print(paste("Control Abd Mean: ", mean(control.abd), " | Control Abd SE: ", sd(control.abd)/sqrt(length(control.abd))))

print(paste("WF Rch Mean: ", mean(wf.rch), " | WF Abd SE: ", sd(wf.rch)/sqrt(length(wf.rch))))
print(paste("WF Abd Mean: ", mean(wf.abd), " | WF Abd SE: ", sd(wf.abd)/sqrt(length(wf.abd))))


mes(mean(wf.rch), mean(control.rch), sd(wf.rch), sd(control.rch), length(wf.rch), length(control.rch))
mes(mean(wf.abd), mean(control.abd), sd(wf.abd), sd(control.abd), length(wf.abd), length(control.abd))

#########################
## Castillo and Wagner ##
#########################

############################################################
## Data in Table 4, sample sizes reported in methodology. ##
############################################################

mes(9, 1.75, 1.82, 0.95, 10, 10) ## richness
mes(145, 50.61, 122, 89.01, 10, 10) ## abundance

mes(8, 2.44, 2.25, 0.598, 10, 10)
mes(29.25, 22.44, 151, 101.96, 10, 10)

mes(8.25, 2.98, 1.25, .43, 6, 6)
mes(25, 7.39, 13.5, 3.77, 6, 6)

######################
## Ghandi et al. ##
######################

######################
## Data in Figure 5 ##
######################

mes(7.5, 19, 8.5, 7, 12, 18)
mes(20, 16, 8, 8, 12, 12)

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

####################
## Moretti et al. ##
####################

###################
## Data in Fig 2 ##
###################

mes(8, 8, 0.55*sqrt(24), 0.5*sqrt(24), 24, 24)
mes(90, 95, 30*sqrt(24), 15*sqrt(24), 24, 24)

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

mes(sum(burnedM[,1]), sum(controlsM[,1]), burnedPSD*sqrt(17), controlsPSD*sqrt(15), 17, 15)

############################
## Samu et al.           ##
############################

#####################
## Data in Table 2 ##
#####################

mes(103.6, 32, 38.89*sqrt(20), 14.58*sqrt(20), 20, 20)
mes(9.8, 8.4, 1.69*sqrt(20), 2.55*sqrt(20), 20, 20)

#####################
## Campbell et al. ##
#####################

#####################
## Data in Table 2 ##
#####################

mes(5, 5, 2*sqrt(12), 1*sqrt(12), 12, 12)
mes(104.3, 118, 25.2*sqrt(12), 52.9*sqrt(12), 12, 12)
mes(2, 3, 1*sqrt(12), 1*sqrt(12), 12, 12)
mes(152.3, 106.7, 73*sqrt(12), 44.2*sqrt(12), 12, 12)

####################
## Lewttow et al. ##
####################

#####################
## Data in Table 1 ##
#####################

mes(41, 37, 5, 3, 20, 20) # 2011 abundance
mes(11, 12, 2, 1, 20, 20) # 2011 richness
mes(38, 31, 1, 5, 20, 20) # 2012 abundance
mes(13, 12, 2, 1, 20, 20) # 2012 richness

#####################
## Lopresti et al. ##
#####################

