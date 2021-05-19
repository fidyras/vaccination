
ModelVacc <- function(paramVac, paramStruc, paramDemo, paramSim, repID){
# Simulation parameter
  tMax = paramSim$ptMax;
  nbSimul = paramSim$nbSimul; # this is obsolete. Argument not used 
  byAge = paramSim$byAge;

  # Population data
  propAge = popStruc$pPropAge;
  n = length(propAge);
  matrixCtctAge = popStruc$pMatrix;
  popSize = popStruc$pPopSize;
  pAi = popStruc$pAi
  pIi = popStruc$pIi

  # Demographic parameters
  transmRate = paramDemo$pTransmRate
  propAsympto = paramDemo$pPropAsympto
  epsilon = paramDemo$pEpsilon
  sigma = paramDemo$pSigma

  severityAge = paramDemo$pSeverityAge;
  lethalityAge = paramDemo$pLethalityAge;
  lethalityAgeUnNoticed = paramDemo$pLethalityAgeUnNoticed;
  # n=length(pPropAge); # of age classes

  # # NPI parameters
  # pisol = paramNPI$pPisol;
  # timeIsol = paramNPI$pTimePisol;
  # confinementBegin = paramNPI$pConfinementBegin;
  # confinementDuration = paramNPI$pConfinementDuration;
  # confinementEfficiency = paramNPI$pConfinementEfficiency;
  # confinementTarget = paramNPI$cpConfinementTarget;

  # vaccination parameters
  vaccinationNumber = paramVac$pVaccinationNumber
  vaccinationPerDay = paramVac$pVaccinationPerDay
  vaccinationBegin = paramVac$pVaccinationBegin
  vaccinationAccept = paramVac$pvaccinationAccept
  vaccinationEfficency = paramVac$pVaccinationEfficency
  # vaccinationPerAge = paramVac$pVaccinationPerAge

  
  tau=0.5;
  lag = 6; # vaccinated started to be immune after 6 days
  
  time=c();
  age=c();
  status=c();
  Number=c();
  Simul=c();
  
  # for(iSimul in c(nbSimul)){
    
    cat('\r',repID, "\n")
    flush.console()

    set.seed(repID)

    tCurrent=0;
    #initial conditions
    S=rep(1,n)*popSize*propAge; # susceptible
    E=rep(0,n); # exposed
    A=pAi; # infectious asymptomatic
    I=pIi; # infectious sympotmatic
    R=rep(0,n); # recovered from I
    U=rep(0,n); # recovered from A
    D=rep(0,n); # death from I
    M=rep(0,n); # serious case from I
    NbCases=rep(0,n); # new case per day
    NVac = rep(0,n) # number of doses per age
    
    onVac = 0 # initialize
    ovac = rep(0,n) #
    
    t=0;
    while(t<tMax){
      #Main body of equations
      S[which(S<0)]=0;
      E[which(E<0)]=0;
      I[which(I<0)]=0;
      A[which(A<0)]=0;
      R[which(R<0)]=0;
      U[which(U<0)]=0;
      M[which(M<0)]=0;
      D[which(D<0)]=0;

      # Do vaccination
      if(vaccinationNumber > 0 & t > (vaccinationBegin + lag)){
        nVac = min(vaccinationPerDay*tau, vaccinationNumber)
        vac = distributeVac(nVac, S, vaccinationAccept)
        # cat(vac, "\n")
        for(i in 1:n){
          temp = rbinom(1, vac[i], vaccinationEfficency[i])
          S[i] = floor(S[i] - temp)
          R[i] = floor(R[i] + temp)
        }
        vaccinationNumber = vaccinationNumber - sum(vac)
        
        onVac <- onVac + nVac # updating as export is based on days not tau step
        ovac <- ovac + vac
      }

      #Calculating force of direct transmission, based on age contact network
      lambda=rep(0,n);
      for(i in 1:n){
        for(j in 1:n){
          beta1=transmRate*matrixCtctAge[j,i];
          # for(k in 1:length(confinementBegin)){
          #   if(t>=confinementBegin[k] & t<(confinementBegin[k]+confinementDuration[k])){
          #     if(j==confinementTarget[k]){
          #       beta1=beta1*(1-confinementEfficiency[k]);
          #     }
          #   }
          # }
          # if(t>=timeIsol){
          #   lambda[i]=lambda[i]+beta1*S[i]*( (1-pisol)*I[j]+(1-pisol)*A[j]);
          # }
          # else{
            lambda[i]=lambda[i]+beta1*S[i]*(I[j]+A[j]);
          # }
        }
      }
      
      
      for(i in 1:n){ #for each age group
        #Rates calculation
        rates=c();
        rates=c(rates,lambda[i]); #S-1,E+1
        rates=c(rates,(1-propAsympto)*epsilon*E[i] ); #E-1,I+1
        rates=c(rates,propAsympto*epsilon*E[i]);
        rates=c(rates,sigma*A[i]); #A-1,U+1
        rates=c(rates,sigma*(1-lethalityAge[i])*(1-severityAge[i])*I[i]); #I-1,R+1
        rates=c(rates,sigma*lethalityAge[i]*I[i]); #I-1,D+1
        rates=c(rates,sigma*severityAge[i]*I[i]); #I-1,M+1
        rates=c(rates,sigma*lethalityAgeUnNoticed[i]*A[i]); #A-1,D+1
        
        #Rates application
        Num=rpois(1,rates[1]*tau);
        S[i]=S[i]-Num;
        E[i]=E[i]+Num;
        Num=rpois(1,rates[2]*tau);
        E[i]=E[i]-Num;
        I[i]=I[i]+Num;
        NbCases[i]=NbCases[i]+Num;
        Num=rpois(1,rates[3]*tau);
        E[i]=E[i]-Num;
        A[i]=A[i]+Num;#E-1,A+1
        Num=rpois(1,rates[4]*tau);
        A[i]=A[i]-Num;
        U[i]=U[i]+Num;
        Num=rpois(1,rates[5]*tau);
        I[i]=I[i]-Num;
        R[i]=R[i]+Num;
        Num=rpois(1,rates[6]*tau);
        I[i]=I[i]-Num;
        D[i]=D[i]+Num;
        Num=rpois(1,rates[7]*tau);
        I[i]=I[i]-Num;
        M[i]=M[i]+Num;
        Num=rpois(1,rates[8]*tau);
        A[i]=A[i]-Num;
        D[i]=D[i]+Num;
      }
      
      #Saving data
      if(t>tCurrent){
        if(byAge==FALSE){
          # number in each status
          time=c(time,rep(tCurrent,8));
          age=c(age,rep(0,8));
          status=c(status,c("S","E", "A","U", "I", "R", "D", "M"));
          Number=c(Number,sum(S),sum(E),sum(A),sum(U),sum(I),sum(R),sum(D),sum(M));
          Simul=c(Simul,rep(repID,8));

          # number of new cases
          time=c(time,tCurrent);
          age=c(age,0);
          status=c(status,"newC");
          Number=c(Number,sum(NbCases));
          Simul=c(Simul,repID);
          NbCases = rep_len(0, n)

          # number of doses
          time=c(time,tCurrent);
          age=c(age,0);
          status=c(status,"nVac");
          Number=c(Number, onVac);
          Simul=c(Simul,repID);
          onVac = 0;
          ovac = rep_len(0,n);
        }else{
          time = c(time, rep(tCurrent, 8 * n ))
          age = c(age, rep(1:n, 8))
          status = c(status, rep(c("S","E", "A","U", "I", "R", "D", "M"), each = n))
          Number= c(Number, c(S, E, A, U, I, R, D, M))
          Simul = c(Simul, rep(repID, 8 * n))

          # number of new cases
          time=c(time,rep(tCurrent, n));
          age=c(age,1:n);
          status=c(status,rep("newC", n));
          Number=c(Number, NbCases); 
          Simul=c(Simul,rep(repID, n));
          NbCases = rep_len(0, n); # re-initialization

          # number of doses
          time=c(time,rep(tCurrent, n));
          age=c(age,1:n);
          status=c(status,rep("nVac", n));
          Number=c(Number, ovac); 
          Simul=c(Simul,rep(repID, n));
          onVac = 0;
          ovac = rep_len(0,n);
        } 
        tCurrent=tCurrent+1;
      }
      t=t+tau;
    } # end loop for one simulation

    time=c(time,rep(tCurrent,n));
    age=c(age,1:n);
    status=c(status,rep("allD",n));
    Number=c(Number,D);
    Simul=c(Simul,rep(repID,n));
    time=c(time,rep(tCurrent,n));
    age=c(age,1:n);
    status=c(status,rep("allR",n));
    Number=c(Number,R);
    Simul=c(Simul,rep(repID,n));
  # } # end loop for replicates

  returnDf=data.frame(time=time,
                      age=age,
                      status=status,
                      Number=Number,
                      Simul=Simul);
  return(returnDf)
}

# compute the Number people vaccinated in each age class
distributeVac <- function(nVac, S, pAccept){
  ll <- length(S)
  result <- rep(0,ll)
  nV <- nVac
  for(i in ll:1){
    if(nV > 0){
      acceptS <- floor(pAccept[i] * S[i])
      if(nV <= acceptS){ # if not enough for the age class
        result[i] <- nV
        nV <- 0
      } else{ # if there are more vaccines available for the class
        result[i] <- acceptS
        nV <- nV - acceptS
      }
    }
  }
  return(result)
}
