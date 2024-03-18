#MCMC for all components
set.seed(1223)
source("Simulation.R")
source("loglikelihood.R")

SimulatedData<- SimulationResults[[1]]
init.density<- c(0.6666667, 0.3333333)
e.it<- 1000
R<- -1*Westmidlands_adjmat
diag(R)<- -rowSums(R, na.rm = T)
qr(R)$rank

#Prior density for Trend components (r_t)
randomwalk2<- function(componentR, PrecisionR){
  time<- length(componentR)
  Sumres<- 0
  for(i in 3:time){
    res<- (componentR[i-2] - (2 * componentR[i-1]) + componentR[i])^2
    Sumres<- Sumres + res
  }
  return((time - 2)/2 * log(PrecisionR) - PrecisionR/2 * Sumres)
}

#Prior density for Seasonal components (s_t)
seasonalComp<- function(x, z){
  time<- length(x)
  Sumres<- 0
  for(i in 12:time){
    res<- (sum(x[(i-11):(i-0)]))^2
    Sumres<- Sumres + res
  }
  return((time - 11)/2 * log(z) - z/2 * Sumres)
}

#Intrinsic GMRF density for spatial components (u_i)
logIGMRF1<- function(x, y, z) {
  n = nrow(z)
  sumC = sum(x[1:(n-1)])
  x = c(x[1:(n-1)], -sumC)
  return ((n - 1)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
}

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

num_iteration<- 150000
MC_chain<- matrix(NA, nrow=num_iteration, ncol=155)
MC_chain[1,]<- c(runif(1), runif(1), kappa_r, kappa_s, kappa_u, r, s, u)
accepted<- 0
likelihoodcurrent<- loglikelihood(SimulatedData,MC_chain[1,6:65],MC_chain[1,66:125],MC_chain[1,126:155],G(MC_chain[1,1],MC_chain[1,2]),init.density,e.it)
zigmaR<- diag(rep(0.0000008, time), nrow = time, ncol = time)
zigmaS<- diag(rep(0.0000008, time), nrow = time, ncol = time)
zigmaU<- diag(rep(0.0001, ndept-1), nrow=ndept-1, ncol=ndept-1)
optconstantRS<- 2.38^2/time
optconstantU<- 2.38^2/(ndept-1)

for(i in 2:num_iteration){
  proposedG12<- abs(rnorm(1,mean=MC_chain[i-1,1], sd=0.01))
  if(proposedG12>1) proposedG12=2-proposedG12
  proposedG21<- abs(rnorm(1,mean=MC_chain[i-1,2], sd=0.01))
  if(proposedG21>1) proposedG21=2-proposedG21
  
  proposedkappaR<- abs(rnorm(1,mean=MC_chain[i-1,3], sd= 3500))
  proposedkappaS<- abs(rnorm(1,mean=MC_chain[i-1,4], sd= 300))
  proposedkappaU<- abs(rnorm(1,mean=MC_chain[i-1,5], sd= 2.7))
  
  proposedRcomps<- rmvnorm(1, mean = MC_chain[i-1, 6:65], sigma = zigmaR)
  
  proposedScomps<- rmvnorm(1, mean = MC_chain[i-1, 66:125], sigma = zigmaS) 
  
  proposedspatcomps<- rmvnorm(1, mean=MC_chain[i-1, 126:154], sigma = zigmaU) 
  proposedspatcomps<- c(proposedspatcomps, -sum(proposedspatcomps))
  
  
  priorcurrent12<- dbeta(MC_chain[i-1,1], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed12<- dbeta(proposedG12, shape1 = 1, shape2 = 1, log=TRUE) 
  priorcurrent21<- dbeta(MC_chain[i-1,2], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed21<- dbeta(proposedG21, shape1 = 1, shape2 = 1, log=TRUE)
  
  priorcurrentkappaR<- dgamma(MC_chain[i-1,3], shape = 1, rate = 0.0001, log=TRUE)
  priorproposedkappaR<- dgamma(proposedkappaR, shape = 1, rate = 0.0001, log=TRUE)
  
  priorcurrentkappaS<- dgamma(MC_chain[i-1,4], shape = 1, rate = 0.0001, log=TRUE)
  priorproposedkappaS<- dgamma(proposedkappaS, shape = 1, rate = 0.0001, log=TRUE)
  
  priorcurrentkappaU<- dgamma(MC_chain[i-1,5], shape = 1, rate = 0.01, log=TRUE)
  priorproposedkappaU<- dgamma(proposedkappaU, shape = 1, rate = 0.01, log=TRUE)
  
  priorcurrentRcomps<- randomwalk2(MC_chain[i-1, 6:65], MC_chain[i-1, 3])
  priorproposedRcomps<- randomwalk2(proposedRcomps, proposedkappaR)
  
  priorcurrentScomps<- seasonalComp(MC_chain[i-1, 66:125], MC_chain[i-1, 4])
  priorproposedScomps<- seasonalComp(proposedScomps, proposedkappaS)
  
  priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 126:155], MC_chain[i-1, 5], R)
  priorproposedUcomps<- logIGMRF1(proposedspatcomps, proposedkappaU, R)
  
  proposalproposedRcomps<- sum(dmvnorm(proposedRcomps, mean = MC_chain[i-1, 6:65], sigma = zigmaR, log = TRUE))
  proposalcurrentRcomps<- sum(dmvnorm(MC_chain[i-1, 6:65], mean = proposedRcomps, sigma = zigmaR, log = TRUE))
  
  proposalproposedScomps<- sum(dmvnorm(proposedScomps, mean = MC_chain[i-1, 66:125], sigma = zigmaS, log = TRUE))
  proposalcurrentScomps<- sum(dmvnorm(MC_chain[i-1, 66:125], mean = proposedScomps, sigma = zigmaS, log = TRUE))
  
  proposalproposedcompsU<- sum(dmvnorm(proposedspatcomps[-30], mean = MC_chain[i-1, 126:154], sigma = zigmaU, log = TRUE))
  proposalcurrentcompsU<- sum(dmvnorm(MC_chain[i-1, 126:154], mean = proposedspatcomps[-30], sigma = zigmaU, log = TRUE))
  
  likelihoodproposed<- loglikelihood(SimulatedData,proposedRcomps,proposedScomps,proposedspatcomps,G(proposedG12,proposedG21),init.density,e.it)
  
  mh.ratio<- exp(likelihoodproposed + priorproposed12  + priorproposed21 + priorproposedkappaR + priorproposedkappaS + priorproposedkappaU + priorproposedRcomps + priorproposedScomps + priorproposedUcomps + proposalcurrentRcomps + proposalcurrentScomps + proposalcurrentcompsU
                  - likelihoodcurrent - priorcurrent12 - priorcurrent21 - priorcurrentkappaR - priorcurrentkappaS - priorcurrentkappaU - priorcurrentRcomps - priorcurrentScomps - priorcurrentUcomps - proposalproposedRcomps - proposalproposedScomps - proposalproposedcompsU)
  
  print(mh.ratio)
  
  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    MC_chain[i,1]<- proposedG12
    MC_chain[i,2]<- proposedG21
    MC_chain[i,3]<- proposedkappaR
    MC_chain[i,4]<- proposedkappaS
    MC_chain[i,5]<- proposedkappaU
    MC_chain[i,6:65]<- proposedRcomps
    MC_chain[i,66:125]<- proposedScomps
    MC_chain[i,126:155]<- proposedspatcomps
    
    likelihoodcurrent<- likelihoodproposed
    
    accepted<- accepted + 1
  }
  else{
    MC_chain[i,1]<- MC_chain[i-1,1]
    MC_chain[i,2]<- MC_chain[i-1,2]
    MC_chain[i,3]<- MC_chain[i-1,3]
    MC_chain[i,4]<- MC_chain[i-1,4]
    MC_chain[i,5]<- MC_chain[i-1,5]
    MC_chain[i,6:65]<- MC_chain[i-1,6:65]
    MC_chain[i,66:125]<- MC_chain[i-1,66:125]
    MC_chain[i,126:155]<- MC_chain[i-1,126:155]
    
    accepted<- accepted + 0
  }
  
  #Adapting zigmaR
  if(i==50){
    epsilonR<- 0.7
    XnR<- MC_chain[1:i, 6:65]
    XnbarR <- colMeans(XnR) 
    zigmaR <- cov(XnR) + epsilonR*diag(rep(1, time))
    zigmaR<- optconstantRS * zigmaR
  } else if (i > 50){ 
    XnbarPrevR <- XnbarR
    XnbarR <- (i*XnbarR + MC_chain[i, 6:65])/(i+1)
    zigmaR <- ((i-1)*zigmaR + tcrossprod(MC_chain[i, 6:65]) + i*tcrossprod(XnbarPrevR) - (i+1)*tcrossprod(XnbarR) + epsilonR*diag(rep(1,time)))/i
    zigmaR<- optconstantRS * zigmaR
    #print(zigmaR)
  }   
  
  #Adapting zigmaS
  if(i==50){
    epsilonS<- 0.7
    XnS<- MC_chain[1:i, 66:125]
    XnbarS <- colMeans(XnS) 
    zigmaS <- cov(XnS) + epsilonS*diag(rep(1, time))
    zigmaS<- optconstantRS * zigmaS
  } else if (i > 50){ 
    XnbarPrevS <- XnbarS
    XnbarS <- (i*XnbarS + MC_chain[i, 66:125])/(i+1)
    zigmaS <- ((i-1)*zigmaS + tcrossprod(MC_chain[i, 66:125]) + i*tcrossprod(XnbarPrevS) - (i+1)*tcrossprod(XnbarS) + epsilonS*diag(rep(1,time)))/i
    zigmaS<- optconstantRS * zigmaS
    #print(zigmaS)
  }   
  
  #Adapting zigmaU
  if(i==50){
    epsilonU<- 0.07
    XnU<- MC_chain[1:i, 126:154]
    XnbarU <- colMeans(XnU) 
    zigmaU <- cov(XnU) + epsilonU*diag(rep(1, ndept-1))
    zigmaU<- optconstantU * zigmaU
  } else if (i > 50){ 
    XnbarPrevU <- XnbarU
    XnbarU <- (i*XnbarU + MC_chain[i, 126:154])/(i+1)
    zigmaU <- ((i-1)*zigmaU + tcrossprod(MC_chain[i, 126:154]) + i*tcrossprod(XnbarPrevU) - (i+1)*tcrossprod(XnbarU) + epsilonU*diag(rep(1,ndept-1)))/i
    zigmaU<- optconstantU * zigmaU
    #print(zigmaU)
  } 
} 

#Results
RsMCMC<- numeric(time)
for(i in 1:time){
  RsMCMC[i] = mean(MC_chain[-(1:20000),5+i])
}
RsMCMC
r

SsMCMC<- numeric(time)
for(i in 1:time){
  SsMCMC[i] = mean(MC_chain[-(1:20000),65+i])
}
SsMCMC
s

UsMCMC<- numeric(ndept)
for(i in 1:ndept){
  UsMCMC[i] = mean(MC_chain[-(1:20000),125+i])
}
UsMCMC
u

colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:60, sep=""), paste("s", 1:60, sep=""), paste("u", 1:30, sep="")))

data_df <- as.data.frame(MC_chain)
TrueValues<- c(0.2, 0.4, kappa_r, kappa_s, kappa_u, r, s, u)

#Histograms
par(mfrow=c(3, 3))  
for (i in 1:ncol(data_df)) {
  hist(data_df[-(1:2000), i], main = colnames(data_df)[i], xlab ="", col = "white", border = "black")
  abline(v=TrueValues[i], col="red", lwd=2,lty=1)
}

#Traceplots
par(mfrow=c(3, 3))  
for (i in 1:ncol(data_df)) {
  plot(data_df[, i], type = "l", main = colnames(data_df)[i], xlab ="MCMC iterations", ylab = "", col = "purple")
  abline(h=TrueValues[i], col="red", lwd=2,lty=1)
}

accepted/num_iteration