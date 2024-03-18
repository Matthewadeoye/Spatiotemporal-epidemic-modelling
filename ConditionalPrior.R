#MCMC for trend components using block updates
set.seed(122)
source("Simulation.R")
source("loglikelihood.R")

SimulatedData<- SimulationResults[[1]]
init.density<- c(0.6666667, 0.3333333)
e.it<- 1000
RW2PrecMat<- matrix(0, nrow=60, ncol=60)
RW2PrecMat[1,(1:3)]<- c(1,-2,1)
RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
RW2PrecMat[59,(57:60)]<- c(1,-4,5,-2)
RW2PrecMat[60,(58:60)]<- c(1,-2,1)
for(i in 3:(time-3)){
  RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
}
RW2PrecMat<- kappa_r * RW2PrecMat

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

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

num_iteration<- 2999
MC_chain<- matrix(NA, nrow=(num_iteration+2), ncol=time)
MC_chain[1,]<- r
ACCEPTED<- 0
likelihoodcurrent<- loglikelihood(SimulatedData,MC_chain[1, 1:60],s,u,G(0.2, 0.4),init.density,e.it)
zigma<- RW2PrecMat
diag(zigma)<- diag(zigma) + 0.000000001
zigma<- solve(zigma)

for(i in seq(from = 2, to = num_iteration, by = 3)){
  for(a in 1:3){
  if(a==1){ 
    conditionalmean<- -solve(RW2PrecMat[1:20, 1:20]) %*% RW2PrecMat[1:20, 21:60] %*% MC_chain[i-1, 21:60]
    conditionalprec<- RW2PrecMat[1:20, 1:20]
    conditionalcov<- solve(conditionalprec)
    proposedRcomps<- rmvnorm(1, mean = conditionalmean, sigma = conditionalcov) 
    proposedRcomps<- c(proposedRcomps, MC_chain[i-1, 21:60])
    
    likelihoodproposed<- loglikelihood(SimulatedData,proposedRcomps,s,u,G(0.2, 0.4),init.density,e.it)
    
    priorcurrentRcomps<- randomwalk2(MC_chain[i-1, 1:60], kappa_r)
    priorproposedRcomps<- randomwalk2(proposedRcomps, kappa_r)
    
    #proposalproposedRcomps<- sum(dmvnorm(proposedRcomps, mean = MC_chain[i-1, 1:60], sigma = zigma, log = TRUE))
    #proposalcurrentRcomps<- sum(dmvnorm(MC_chain[i-1, 1:60], mean = proposedRcomps, sigma = zigma, log = TRUE))
    
    mh.ratio<- exp(likelihoodproposed + priorproposedRcomps - likelihoodcurrent - priorcurrentRcomps)
    print(mh.ratio)
    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      ACCEPTED<- ACCEPTED + 1
      MC_chain[i,1:60]<- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
    }
    else{
      MC_chain[i, 1:60]<- MC_chain[i-1, 1:60]
      ACCEPTED<- ACCEPTED + 0
    }
  }
    
  else if(a==2){ 
    conditionalmean<- -solve(RW2PrecMat[21:40, 21:40]) %*% (RW2PrecMat[21:40, 1:20]%*%MC_chain[(i+1)-1, 1:20] + RW2PrecMat[21:40, 41:60]%*%MC_chain[(i+1)-1, 41:60]) 
    conditionalprec<- RW2PrecMat[21:40, 21:40]
    conditionalcov<- solve(conditionalprec)
    proposedRcomps<- rmvnorm(1, mean = conditionalmean, sigma = conditionalcov) 
    proposedRcomps<- c(MC_chain[(i+1)-1, 1:20], conditionalmean, MC_chain[(i+1)-1, 41:60]) 
    
    likelihoodproposed<- loglikelihood(SimulatedData,proposedRcomps,s,u,G(0.2, 0.4),init.density,e.it)
    
    priorcurrentRcomps<- randomwalk2(MC_chain[(i+1)-1, 1:60], kappa_r)
    priorproposedRcomps<- randomwalk2(proposedRcomps, kappa_r)
    
    #proposalproposedRcomps<- sum(dmvnorm(proposedRcomps, mean = MC_chain[(i+1)-1, 1:60], sigma = zigma, log = TRUE))
    #proposalcurrentRcomps<- sum(dmvnorm(MC_chain[(i+1)-1, 1:60], mean = proposedRcomps, sigma = zigma, log = TRUE))
    
    
    mh.ratio<- exp(likelihoodproposed + priorproposedRcomps - likelihoodcurrent - priorcurrentRcomps)
    print(mh.ratio)
    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      ACCEPTED<- ACCEPTED + 1
      MC_chain[(i+1),1:60]<- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
    }
    else{
      MC_chain[(i+1), 1:60]<- MC_chain[(i+1)-1, 1:60]
      ACCEPTED<- ACCEPTED + 0
    }
  }

  else if(a==3){ 
    conditionalmean<- -solve(RW2PrecMat[41:60, 41:60]) %*% RW2PrecMat[41:60, 1:40] %*% MC_chain[(i+2)-1, 1:40]
    conditionalprec<- RW2PrecMat[41:60, 41:60]
    conditionalcov<- solve(conditionalprec)
    proposedRcomps<- rmvnorm(1, mean = conditionalmean, sigma = conditionalcov)
    proposedRcomps<- c(MC_chain[(i+2)-1, 1:40], conditionalmean) 
    
    likelihoodproposed<- loglikelihood(SimulatedData,proposedRcomps,s,u,G(0.2, 0.4),init.density,e.it)
    
    priorcurrentRcomps<- randomwalk2(MC_chain[(i+2)-1, 1:60], kappa_r)
    priorproposedRcomps<- randomwalk2(proposedRcomps, kappa_r)
    
    #proposalproposedRcomps<- sum(dmvnorm(proposedRcomps, mean = MC_chain[(i+2)-1, 1:60], sigma = zigma, log = TRUE))
    #proposalcurrentRcomps<- sum(dmvnorm(MC_chain[(i+2)-1, 1:60], mean = proposedRcomps, sigma = zigma, log = TRUE))
    
    
    mh.ratio<- exp(likelihoodproposed + priorproposedRcomps - likelihoodcurrent - priorcurrentRcomps)
    print(mh.ratio)
    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      ACCEPTED<- ACCEPTED + 1
      MC_chain[(i+2),1:60]<- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
    }
    else{
      MC_chain[(i+2), 1:60]<- MC_chain[(i+2)-1, 1:60]
      ACCEPTED<- ACCEPTED + 0
    }
   } 
  }
}

RsMCMC<- numeric(time)
for(i in 1:time){
  RsMCMC[i] = mean(MC_chain[-((num_iteration+1):(num_iteration+3)), i])
}
RsMCMC
r

colnames(MC_chain) <- paste("r", 1:60, sep="")

data_df <- as.data.frame(MC_chain)
TrueValues<- r

#Histograms
par(mfrow=c(3, 3))  
for (i in 1:ncol(data_df)) {
  hist(data_df[-(1:20), i], main = colnames(data_df)[i], xlab ="", col = "white", border = "black")
  abline(v=TrueValues[i], col="red", lwd=2,lty=1)
}

#Traceplots
par(mfrow=c(3, 3))  
for (i in 1:ncol(data_df)) {
  plot(data_df[, i], type = "l", main = colnames(data_df)[i], xlab ="MCMC iterations", ylab = "", col = "purple")
  abline(h=TrueValues[i], col="red", lwd=2,lty=1)
}

