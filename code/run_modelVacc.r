source("code/ModelVacc.R")

paramDATA <- read.csv("data/paramDemo.csv")
pLA <- paramDATA$pLethalityAge
pSA <- paramDATA$pSeverityAge
pLU <- rep(0,7) #baseline scenario: Asymptomatics do not die of COVID #paramDATA$pLethalityAgeUnNoticed 
baseTR <- 1.21e-9 # baseline transmission rate used in Evans

pMatrix <- read.csv("data/ContactMatrix.csv")

tMax <- 365 #one year
tau <- 0.5 #timesteps
nbSimul <- 1
nn <- 7 # number of age class
mm<-nn*tMax*10 #(10:number of compartments SEAIRUMD+nVac+newC)
# basic parameter for simulation
paramSim <- list(ptMax = tMax, nbSimul = nbSimul, byAge = T)

# demographic parameter
paramDemo <- list(pTransmRate = baseTR, pPropAsympto = 0.4, pEpsilon = 1/3, pSigma = 1/5, pLethalityAge = pLA, pSeverityAge = pSA, pLethalityAgeUnNoticed = pLU)

# # parameter for npi (if needed)
# paramNPI <- list(pPisol = 0, pTimePisol = Inf, pConfinementBegin = rep(Inf, nn),pConfinementDuration = rep(Inf, nn), pConfinementEfficiency = rep(Inf, nn), cpConfinementTarget = 1:nn)

# DEFINE REMAINING PARAMETERS 

POPDATA <- read.csv("data/age_dist.csv")
VACDATA <- read.csv("data/vacc_allocation.csv")
REG <- unique(POPDATA$key)
TOTALPOP <- sum(POPDATA$sum_pop)

TOTALVAC <- c(5000000)# ceiling(c(TOTALPOP*0.01)) # total number of vaccines #baseline scenario 25% of TOTALPOPULATION
STARTV <- c(10) # when to start vaccinating
VPERDAY <- ceiling(c(1*20))# proportion of staff * how many they can do in a day
VACCOPT <- 1:1 #c("freq_pop", "freq_60", "freqcases", "freqdeaths","uniform")
VA <-  c(0,0.01,0.1,0.5,1) # acceptance #baseline 0.70 (Transparency international)
VE <- c(1) # vaccine efficiency #baseline scenario: 0.9
RR <- 3:3 # 1:22
REPID <- 1:10

summaryD <- c()
ts<-c()


## set pPropR: initional the proportion of recovered from previous infection
pPropR <- rep(0, nn) 

# RUNNING WITH ALL COMBINATIONS
for(tVac in TOTALVAC){ # total number of vaccines
	for(startVac in STARTV){ # beginning of vaccination 
		for(vpd in VPERDAY){ # proportion of staff willing to do the vaccine and how many per day
			for(vo in VACCOPT){ # allocate vaccine according to size, old, case or death
				for(va in VA){ # vaccine acceptance
					pvA <- rep(va, nn) # assume does not vary with age for now
					for(ve in VE){ # vaccine efficacy 
						pvEff <- rep(ve, nn) # efficacy does not vary with age for now
						for(r in RR){

							# construct population data
							temp <- POPDATA[which(POPDATA$key == REG[r]),]
							pPopSize <- sum(temp$sum_pop) # population size in the region
							pPropAge <- temp$sum_pop/pPopSize #  convert to proportion per age class
							paramDemo$pTransmRate = baseTR * TOTALPOP/pPopSize # scale so that R0 ~ 2.5

							#start with appropriate numbers of infected in each age class based on age distribution
							#start with ten casesvap
							pIi <- round(10 * pPropAge,0)

							#number of asymptomatics
							pAi <- round(pIi/(1 - paramDemo$pPropAsympto),0)

							# parameters for the structure of the population
							popStruc <- list(pMatrix = pMatrix, pPropAge = pPropAge, pPopSize = pPopSize, pIi = pIi, pAi = pAi, pPropR = pPropR)

							# construct vaccination parameter
							pp <- VACDATA[which(VACDATA$key == REG[r]), c(2 + vo,8)] # extract allocation option and the number of medical staff
							nVac <- round(tVac * pp[1], 0)
							vacPerDay <- round(pp[2] * vpd, 0)

							# parameters for vaccination
							paramVac <- list(pVaccinationNumber = nVac, pVaccinationPerDay = vacPerDay, pVaccinationBegin = startVac, pvaccinationAccept = pvA, pVaccinationEfficency = pvEff)
							
							# run the model
							for(repID in REPID){
								output <- ModelVacc(paramVac, paramStruc, paramDemo, paramSim, repID)

								# add the parameter to the output and append the region
								summaryD0 <- tail(output,14)[1:nn, 2:5]
								timeseries <- output[1:mm, 1:5]
								# print(summaryD0)
								regID <- rep(REG[r],nn)
								regNV <- as.numeric(rep(nVac,nn))
								startV <- rep(startVac,nn)
								vPerD <- as.numeric(rep(vacPerDay,nn))
								allocV <- rep(vo,nn)
								acceptV <- pvA
								totalvac<- tVac
								# cat(pvA, "\n")
								effV <- pvEff
								
								# output of the timeseries for all regions
								regID_ts <- rep(REG[r],mm)
								allocV_ts <- rep(vo,mm)
								timeseries<-cbind(timeseries,regID_ts,allocV_ts,acceptV)
								 ts<-rbind(timeseries,ts)
								
								summaryD0 <- cbind(summaryD0, regID, regNV, startV, vPerD, vpd, vo, acceptV, effV,totalvac) 
								summaryD <- rbind(summaryD, summaryD0)

								outname <- paste("Simdata/sim_vacc_reg",REG[r],"vAlloc", vo, "tv", tVac, "dv", vpd, "sv", startVac,"ev", ve, "av", va , "rep", repID, ".csv", sep = "_")
								write.csv(output, outname)
							}
						}
					}	
				}
			}
		}
	}
}


outname_ts <- paste("Simdata/SEIR","accept",va, "start", startVac, ".csv", sep = "_")
write.csv(ts,outname_ts)
outname <- paste("Simdata/visual.csv")
write.csv(summaryD, outname)
