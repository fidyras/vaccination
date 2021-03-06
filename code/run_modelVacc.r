source("code/ModelVacc.R")

paramDATA <- read.csv("data/paramDemo.csv")
pLA <- 10*paramDATA$pLethalityAge
pSA <- paramDATA$pSeverityAge
pLU <- rep(0,7) #baseline scenario: Asymptomatics do not die of COVID #paramDATA$pLethalityAgeUnNoticed 
baseTR <- 1.51e-9 # baseline transmission rate used in Evans

pMatrix <- read.csv("data/ContactMatrix.csv")

tMax <- 365 #365 #one year
tau <- 0.5 #timesteps
nbSimul <- 1
nn <- 7 # number of age class
mm<-nn*tMax*10 #(10:number of compartments SEAIRUMD+nVac+newC)
# basic parameter for simulation
paramSim <- list(ptMax = tMax, nbSimul = nbSimul, byAge = T)

# waning immunity rate: the same for all state 
pW <- 1/182 # roughly half a year. 

# demographic parameter
paramDemo <- list(pTransmRate = baseTR, pPropAsympto = 0.4, pEpsilon = 1/3, pSigma = 1/5, pLethalityAge = pLA, pSeverityAge = pSA, pLethalityAgeUnNoticed = pLU, pWane = pW)

# # parameter for npi (if needed)
# paramNPI <- list(pPisol = 0, pTimePisol = Inf, pConfinementBegin = rep(Inf, nn),pConfinementDuration = rep(Inf, nn), pConfinementEfficiency = rep(Inf, nn), cpConfinementTarget = 1:nn)

# DEFINE REMAINING PARAMETERS 

POPDATA <- read.csv("data/age_dist.csv")
VACDATA <- read.csv("data/vacc_allocation.csv")
REG <- unique(POPDATA$key)
TOTALPOP <- sum(POPDATA$sum_pop)

TOTALVAC <-  ceiling(c(TOTALPOP*0.2)) # total number of vaccines #baseline scenario 25% of TOTALPOPULATION
STARTV <- c(10) # when to start vaccinating
VPERDAY <- ceiling(c(0.5*20))# proportion of staff * how many they can do in a day
VACCOPT <- 1:1 #c("freq_pop", "freq_60", "freqcases", "freqdeaths","uniform")
VA <-  c(0,0.7) # acceptance #baseline 0.70 (Transparency international)
VE <- c(0.76) # vaccine efficiency #baseline scenario: 0.9
RR <- 1:22 # 1:22
REPID <- 1:25
r0fact<-1 #factor by which to increase/decrease the r0 (initial is 2.5)
 r0<-r0fact*2.5

testing <- F # test for presence of antibodies on site (vaccinate only seronegative)
vacByAge <- T # whether to prioritize older class or just distribute randomly 
onlyS <- T # only susceptible come to be vaccinated
wscale <- 0 # change to 1 for 6 months and 2 for 3 months

#######set baseline seroprevalence by age and region
spG<-c(0) #national seroprevalence
## set pPropR: initional the proportion of recovered from previous infection
NR <-  matrix(0,nrow=22,ncol=nn) 
regR<-round(spG*TOTALPOP*VACDATA$freq_pop,0)# number of seropositives per region distribution
spAge<-c(0.05,0.1,0.2,0.2,0.2,0.15,0.1) #change seroprevalence by age class
for (rr in 1:22){
  NR[rr,] <- round(regR[rr]*spAge)
}

summaryD <- c()
ts<-c()

# RUNNING WITH ALL COMBINATIONS
system.time(
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
							paramDemo$pTransmRate = r0fact*baseTR * TOTALPOP/pPopSize # scale so that R0 ~ 2.5
							paramDem$pWane = wscale*paramDem$pWane 

							#start with appropriate numbers of infected in each age class based on age distribution
							#either start with ten cases or 1 infected individual per 100 000 people in the region
 							pIi <- round(round(sum(temp$sum_pop/100000)) * pPropAge,0)

							#number of asymptomatics
							pAi <- round(pIi/(1 - paramDemo$pPropAsympto),0)

							# parameters for the structure of the population
							popStruc <- list(pMatrix = pMatrix, pPropAge = pPropAge, pPopSize = pPopSize, pIi = pIi, pAi = pAi, pnR = NR[r,])

							# construct vaccination parameter
							pp <- VACDATA[which(VACDATA$key == REG[r]), c(2 + vo,8)] # extract allocation option and the number of medical staff
							nVac <- round(tVac * pp[1], 0)
							vacPerDay <- round(pp[2] * vpd, 0)

							# parameters for vaccination
							paramVac <- list(pVaccinationNumber = nVac, pVaccinationPerDay = vacPerDay, pVaccinationBegin = startVac, pvaccinationAccept = pvA, pVaccinationEfficency = pvEff)
							
							# run the model
							for(repID in REPID){
								output <- ModelVacc(paramVac, paramStruc, paramDemo, paramSim, testing, vacByAge, onlyS, repID)

								# add the parameter to the output and append the region
								summaryD0 <- tail(output,14)[1:nn, 2:5]
								 timeseries <- output[1:mm, 1:5]
								 #print(summaryD0)
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
								
								summaryD0 <- cbind(summaryD0, regID, regNV, startV, vPerD, vpd, vo, acceptV, effV,totalvac,spG, testing, vacByAge, onlyS, r0) 
								summaryD <- rbind(summaryD, summaryD0)

								#outname <- paste("Simdata/sim_vacc_reg",REG[r],"vAlloc", vo, "tv", tVac, "dv", vpd, "sv", startVac,"ev", ve, "av", va , "rep", repID, ".csv", sep = "_")
								#write.csv(output, outname)
							}
						}
					}	
				}
			}
		}
	}
}
)
	
#outname_ts <- paste("Simdata/SEIR","accept",va, "start", startVac, ".csv", sep = "_")
#write.csv(ts,outname_ts)
#outname <- paste("Simdata/baseline_seroprevalence0.2_uniform.csv")
#write.csv(summaryD, outname)

basedata_r2<-basedata_r%>%
   rbind(summaryD)