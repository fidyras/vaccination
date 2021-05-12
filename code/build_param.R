# BASIC PARAMETERS CONSTRUCTION 


#### Base Model Parameters ####

# DEMOGRAPHIC PARAMETERS per age class
nn <- 7

pLethalityAge <- c(0.000185, 0.00149, 0.0138, 0.0408, 0.0764, 0.344, mean(1.11, 2.45, 3.8))/100 # from ifr_list wuhan
pLethalityAgeUnNoticed <- rep(0.008,nn); #deaths from Asymptomatic class
#this assumes that all those who died were also counted in the severe cases, so we subtract so we don't double count

#rates from Ferguson
hosp.rates <- c(0.1, 0.3, 1.2, 3.2,4.9,10.2, mean(c(15.6, 24.3, 27.3)))/100
#prop of hospitalized patients who need critical care
crit.care.rates <- c(5,5,5,5,6.3,12.2, mean(c(27.4, 43.2,70.9)))/100
severe.rates <- hosp.rates *crit.care.rates
pSeverityAge <- severe.rates - pLethalityAge

# POPULATION STRUCTURE

# contact matrix

# Contact Matrix
# original from Evans
#taken directly from Roche code
ctc <- read.csv("data/mozambique-ctc.csv", header = F);
pMatrix0 <- matrix(0,8,8);
#combine ages from 16 to 8 groups
for(i in 1:8){
  for(j in 1:8){
    pMatrix0[i,j]=sum(ctc[c( (i-1)*2+1,(i-1)*2+2),c( (j-1)*2+1,(j-1)*2+2)]);
  }
}

# Combine so age class 60 and 70 to 60+
pMatrix <- matrix(0,7,7)
pMatrix <- pMatrix0
pMatrix[7,] <- pMatrix[7,] + pMatrix[8,]
pMatrix[,7] <- pMatrix[,7] + pMatrix[,8]
pMatrix <- pMatrix[-8,-8]
 
# epsilon = 1/3; #E to I
# sigma = 1/5; #I to something? *1/infectious period

# pPropAsympto = 0.8;
# pTransmRate = 1.21e-9
# popAge[7] <- sum(popAge[c(7,8)])
# sum(pTransmRate*(popAge[-8]*rowSums(pMatrix)))/(epsilon+sigma) #R0

write.csv(pMatrix, "data/ContactMatrix.csv", row.names = F, col.names = F)

#Initial Conditions
distribAge=pPropAge
#start with appropriate numbers of infected in each age class based on age distribution
#start with ten cases
pIi <- round(10 * distribAge,0)

#number of asymptomatics
pAi <- round(pIi/(1-pPropAsympto),0)

ageClass <- c("0-9","10-19","20-29", "30-39", "40-49", "50-59", "60+")

# EXPORT FIXED PARAMETERS
PARAM <- data.frame(
	ageClass = ageClass, 
	pLethalityAge = pLethalityAge, 
	pLethalityAgeUnNoticed = pLethalityAgeUnNoticed,
	pSeverityAge = pSeverityAge,)

write.csv(PARAM, "data/paramDemo.csv")

##########################
####### TEST MODEL #######
##########################

setwd("/Users/tanjona/Box Sync/projects/corona_mada/Vaccination")

source("code/ModelVacc.R")

paramDATA <- read.csv("data/paramDemo.csv")
pLA <- paramDATA$pLethalityAge
pSA <- paramDATA$pSeverityAge
pLU <- paramDATA$pLethalityAgeUnNoticed

pMatrix <- read.csv("data/ContactMatrix.csv")

tMax <- 365 #one year
tau <- 0.5 #timesteps
nbSimul <- 1
nn <- 7 # number of age class

# basic parameter for simulation
paramSim <- list(ptMax = tMax, nbSimul = nbSimul, byAge = F)

# demographic parameter
paramDemo <- list(pTransmRate = 1.21e-9, pPropAsympto = 0.8, pEpsilon = 1/3, pSigma = 1/5, pLethalityAge = pLA, pSeverityAge = pSA, pLethalityAgeUnNoticed = pLU)

# # parameter for npi (if needed)
# paramNPI <- list(pPisol = 0, pTimePisol = Inf, pConfinementBegin = rep(Inf, nn),pConfinementDuration = rep(Inf, nn), pConfinementEfficiency = rep(Inf, nn), cpConfinementTarget = 1:nn)

# DEFINE REMAINING PARAMETERS 

POPDATA <- read.csv("data/age_dist.csv")
VACDATA <- read.csv("data/vacc_allocation.csv")
REG <- unique(POPDATA$key)
TOTALPOP <- sum(POPDATA$sum_pop)

TOTALVAC <- c(4*10^6) # total number of vaccines
STARTV < c(10) # when to start vaccinating
VPERDAY <- c(0.5*20) # proportion of staff times how many a person can do in a day
VACCOPT <- 1:4 #c("freq_pop", "freq_60", "freqcases", "freqdeaths")
VA <- c(1) # acceptance 
VE <- c(1) # vaccine efficiency


tVac <- 4*10^6
startVac <- 20
vpd <- 0.5*200
vo <- 1
va <- 1; pvA <- rep(va, nn)
ve <- 1; pvEff <- rep(ve, nn)

for(r in 1:1){
	r <- 3
	# construct population data
	temp <- POPDATA[which(POPDATA$key == REG[r]),]
	pPopSize <- sum(temp$sum_pop) # population size in the region
	pPropAge <- temp$sum_pop/pPopSize #  convert to proportion per age class
	paramDemo$pTransmRate = 1.21e-9 * TOTALPOP/pPopSize # scale so that R0 ~ 2.5

	#start with appropriate numbers of infected in each age class based on age distribution
	#start with ten cases
	pIi <- round(10 * pPropAge,0)

	#number of asymptomatics
	pAi <- round(pIi/(1 - paramDemo$pPropAsympto),0)

	# parameters for the structure of the population
	popStruc <- list(pMatrix = pMatrix, pPropAge = pPropAge, pPopSize = pPopSize, pIi = pIi, pAi = pAi)

	# construct vaccination parameter
	pp <- as.numeric(VACDATA[which(VACDATA$key == REG[r]), c(2 + vo, 7)]) # extract allocation option and the number of medical staff
	nVac <- round(tVac * pp[1], 0)
	vacPerDay <- round(pp[2] * vpd, 0)

	# parameters for vaccination
	paramVac <- list(pVaccinationNumber = nVac, pVaccinationPerDay = vacPerDay, pVaccinationBegin = startVac, pvaccinationAccept = pvA, pVaccinationEfficency = pvEff)
	
	# run the model
	output <- ModelVacc(paramVac, paramStruc, paramDemo, paramSim)

	# # add the parameter to the output and append the region
	# summaryD0 <- tail(output,14)[1:nn, 2:5]
	# regID <- rep(REG[r],nn)
	# regNV <- rep(nVac,nn)
	# startV <- rep(startVac,nn)
	# vPerD <- rep(vacPerDay,nn)
	# allocV <- rep(vo,nn)
	# acceptV <- rep(pvA,nn)
	# effV <- pvEff


	# summaryD0 <- cbind(summaryD0, regID, regNV, startV, vPerD, allocV, acceptV, effV) 
	# if(r == 1){
	# 	summaryD <- summaryD0
	# } else{
	# 	summaryD <- rbind(summaryD, summaryD0)
	# } 

	outname <- paste("data/sim_data/test_sim_vacc_reg_",REG[r],".csv", sep = "")
	write.csv(output, outname, row.names = F)
}

# plot results 

data <- output[which(output$status == "S"), c(1,4)]

plot(data$time, data$Number, typ = "l")

tail(output,14)[1:nn, 2:5]