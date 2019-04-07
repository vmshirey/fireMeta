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
