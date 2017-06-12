################################################################################################################################################################################################
#	https://github.com/baktoft/yaps_sciRep
#	Baktoft, Gjelland, Økland and Thygesen (2017). Time-of-arrival based positioning of aquatic animals using YAPS (Yet Another Positioning System).
#	Example code to reproduce YAPS-results presented in figure 5 from the manuscript. Position estimation is based on raw time-of-arrival data obtained from a tow-track conducted in the field.
#	TMB needs to be downloaded and installed beforehand - see https://github.com/kaskr/adcomp/wiki/Download
################################################################################################################################################################################################
#	This example file should be able to run trouble free, once TMB is installed. Tested on Windows 7, using R x64 3.3.1. 
#	In case of problems, please contact: hba@aqua.dtu.dk or submit as an issue on github. Thanks!
################################################################################################################################################################################################
rm(list=ls())
set.seed(42) #for reproducible results

#load required packages - install if necessary
library(zoo)
library(TMB) 
# To test TMB installation
runExample(all=TRUE)

#source support functions
source('supportFuncs.R')

#Position of hydrophones. H14 malfunctioned and is not used - (x,y ) = (-9999,-9999)
Hx = c(0,143.34,105.93,83.45,-86.3,-133.42,-107.29,15.05,36.91,-91.02,-44.82,61.99,113.24, -9999, 86.9,-125.92,21.67,52.63,55.33,-2.12,-41.71,-22.99,-78.85,-83.19,74.96,27.66,44.05,-115.85)
Hy = c(0,-78.85,-85.87,-121.16,22.39,-55.03,31.69,-24.81,-34.47,4.96,6.59,-15.86,-53.76,-9999, -40.88,-33.9,-70.53,-93.06,-46.82,-24.19,-27.04,-60.77,-58.98,-19.39,-66.14,-5.73,-10.58,10.58)

#read and prepare toa-file
toa <- read.table('toa.txt') #Time-of-arrival obtained from the hydrophones.
toa <- t(as.matrix(toa))

#We suggest using a sub-sample for first testing (e.g. first 500 pings) as computation time is faster
#first toa-observation to include
nstart <- 1	
#number of pings to include - increase to 3244 for complete track. Don't go below ~100 as number of detecting hydrophones is quite low in the beginning of the track (often <3 hydrophones detecting each ping), which can cause numerical problems. 
n <- 1000
nmax <- min(n + nstart, ncol(toa))	#last toa-observation to include

toa <- toa[,nstart:nmax]			
T0 <- min(toa, na.rm=TRUE)
toa <- toa - T0	


#Increase pNA (range 0.00 - 1.00 ) to randomly increase probability of missed detections
pNA <- 0.00
toa[which(rbinom(length(toa),1,pNA) == 1)] <- NA

#Subsample the hydrophone array by excluding hydros specified in noHs (range 1 - 28).
noHs <- c(14) #Hydrophone 14 malfunctioned and is not used

#This configuration leaves eight hydrophones to cover the entire study site
# noHs <- c(14,28,16,10,5,8,26,27,19,18,13,3,20,9,25,24,21,22, 15,1)
toa[noHs, ] <- NA


toa[is.na(toa)] <- -9999 			#NA's are translated to -9999

#Get data prepared for TMB
datTmb <- getDatTmb(Hx, Hy, toa)

#Get list of parameters to be estimated by TMB-model
params <- getParams(datTmb)

#Get list of initial values for estimated fixed parameters
inits <- getInits()

#Compile TMB-model - only needed once
compile("yaps.cpp")

#Run YAPS-model and get results - can take a while depending on size of toa
yapsRes <- runYAPS(inits, params, silent=TRUE)

#Plot results and compare to gps and umap

plotRes(yapsRes)

