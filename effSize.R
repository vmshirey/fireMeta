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

####################
## Moretti et al. ##
####################

###################
## Data in Fig 2 ##
###################

mes(8, 8, 0.55*sqrt(24), 0.5*sqrt(24), 24, 24) # richness
mes(90, 95, 30*sqrt(24), 15*sqrt(24), 24, 24) # abundance

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

########################
## Martikainen et al. ##
########################

# values are inverse effects 

mes(103, 88, 19*sqrt(20), 34*sqrt(20), 20, 20) # abundance
mes(5.6, 7.8, 0.9*sqrt(20), 2.2*sqrt(20), 20, 20) # richness

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
## Gongalsky _____ ##
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

b.abun <- c(31, 62, 3, 4, 50, 102, 68, 41, 35)
u.abun <- c(28, 72, 89, 44, 100, 51, 28, 44, 43)

b.rich <- c(7, 10, 3, 3, 5, 10, 9, 7, 5)
u.rich <- c(5, 5, 8, 5, 8, 6, 9, 6, 8)

mes(mean(b.abun), mean(u.abun), sd(b.abun), sd(u.abun), length(b.abun), length(u.abun))
mes(mean(b.rich), mean(u.rich), sd(b.rich), sd(u.rich), length(b.rich), length(u.rich))

## BEES ##

#####################
## Campbell et al. ##
#####################

#####################
## Data in Table 2 ##
#####################

bAbun <- c(104.3, 4.7, 5.0)
bAbunVar <- c(25.2, 2.7, 2.0)
bAbunSD <- sqrt(bAbunVar)
uAbun <- c(118, 4.0, 5.0)
uAbunVar <- c(52.9, 2.0, 1.0)
uAbunSD <- sqrt(uAbunVar)

mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10)

bAbun <- c(152.3, 5.7, 2.0)
bAbunVar <- c(73, 0.88, 1.0)
bAbunSD <- sqrt(bAbunVar)
uAbun <- c(106.7, 17.7, 3.0)
uAbunVar <- c(44.2, 11.2, 1.0)
uAbunSD <- sqrt(uAbunVar)

mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10)

####################
## Lettow et al.  ##
####################

#####################
## Data in Table 1 ##
#####################

mes(41, 37, 5, 3, 20, 20) # 2011 abundance
mes(11, 12, 2, 1, 20, 20) # 2011 richness
mes(38, 31, 1, 5, 20, 20) # 2012 abundance
mes(13, 12, 2, 1, 20, 20) # 2012 richness

#################
## Love et al. ##
#################

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

mes(18, 17, 2.4, 2.0, 3, 3) # richness
mes(469, 161, 70, 32, 3, 3) # abundance

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

########################
## Data in Appendix 1 ##
########################

hi.burn <- c(472, 1314, 818, 607) # abundance
mx.burn <- c(598, 1235, 777, 334)
control <- c(75)

mes(mean(hi.burn), mean(control), sd(hi.burn), 0.0001, 4, 1)
mes(mean(mx.burn), mean(control), sd(mx.burn), 0.0001, 4, 1)

hi.burn <- c(106, 83, 94, 90) # richness
mx.burn <- c(98, 109, 103, 66)
control <- c(25)

mes(mean(hi.burn), mean(control), sd(hi.burn), 0.0001, 4, 1)
mes(mean(mx.burn), mean(control), sd(mx.burn), 0.0001, 4, 1)

####################
## Granath et al. ##
####################

###################
## Data on Dryad ##
###################

control <- c(0,0,0,0,2,2,6,0,9,1,2,2) # abundance
burn <- c(7,7,3,2,2,11,6,2,5,6,4,3)

mes(mean(burn), mean(control), sd(burn), sd(control), 12, 12)

control <- c(0,0,0,0,2,2,2,0,4,1,3,2) # richness
burn <- c(7,6,3,1,2,6,4,2,4,3,2,2)

mes(mean(burn), mean(control), sd(burn), sd(control), 12, 12)

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

## BUTTERFLIES ##

#########################
## Campbell et al 2007 ##
#########################

## Data in table 2

bAbun <- c(13.7, 3.5)
bAbunVar <- c(3.5, 0.5)
bAbunSD <- sqrt(bAbunVar)
uAbun <- c(11.3, 5.3)
uAbunVar <- c(4.9, 1.5)
uAbunSD <- sqrt(uAbunVar)

mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10)

bAbun <- c(20.7, 4.0)
bAbunVar <- c(10.7, 2.5)
bAbunSD <- sqrt(bAbunVar)
uAbun <- c(11.3, 4.0)
uAbunVar <- c(2.7, 0)
uAbunSD <- sqrt(uAbunVar)

mes(mean(bAbun), mean(uAbun), (sum(bAbunSD^2)/10), sum(uAbunSD^2)/10, 10, 10)

#####################
## Fleishman et al ##
#####################

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

mes(2.79, )

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

