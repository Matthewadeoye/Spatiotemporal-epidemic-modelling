#MCMC for seasonal components (s_t)
set.seed(122)
source("Simulation.R")
source("loglikelihood.R")

SimulatedData<- SimulationResults[[1]]
init.density<- c(0.6666667, 0.3333333)
e.it<- 1000

#Prior density for Seasonal components (s_t)
seasonalComp<- function(x, z){
  time<- length(x)
  Sumres<- 0
  for(i in 12:time){
    res<- (sum(x[(i-11):(i-0)]))^2
    Sumres<- Sumres + res
  }
  return((time - 11)/2 * (log(z) - log(2 * pi)) -z/2 * Sumres)
}

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

num_iteration<- 15000
MC_chain<- matrix(NA, nrow=(num_iteration), ncol=time)
MC_chain[1,]<- rep(0, time)
ACCEPTED<- 0
likelihoodcurrent<- loglikelihood(SimulatedData,r,MC_chain[1, 1:60],u,G(0.2, 0.4),init.density,e.it)

zigma<- diag(rep(0.0000008, time), nrow = time, ncol = time)

for(i in 2:num_iteration){
  proposedScomps<- rmvnorm(1, mean = MC_chain[i-1, ], sigma = zigma) 
  
  likelihoodproposed<- loglikelihood(SimulatedData,r,proposedScomps,u,G(0.2, 0.4),init.density,e.it)
  
  priorcurrentScomps<- seasonalComp(MC_chain[i-1, ], kappa_s)
  priorproposedScomps<- seasonalComp(proposedScomps, kappa_s)
  
  proposalproposedScomps<- sum(dmvnorm(proposedScomps, mean = MC_chain[i-1, ], sigma = zigma, log = TRUE))
  proposalcurrentScomps<- sum(dmvnorm(MC_chain[i-1, ], mean = proposedScomps, sigma = zigma, log = TRUE))
  
  mh.ratio<- exp(likelihoodproposed + priorproposedScomps + proposalcurrentScomps - likelihoodcurrent - priorcurrentScomps - proposalproposedScomps)
  print(mh.ratio)
  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    ACCEPTED<- ACCEPTED + 1
    MC_chain[i,1:60]<- proposedScomps
    likelihoodcurrent<- likelihoodproposed
  }
  else{
    MC_chain[i, 1:60]<- MC_chain[i-1, 1:60]
    ACCEPTED<- ACCEPTED + 0
  }
  #Adapting zigma
  if(i==500){
    epsilon<- 0.06
    Xn<- MC_chain[1:i, 1:60]
    Xnbar <- colMeans(Xn) 
    zigma <- cov(Xn) + epsilon*diag(rep(1, time))
    zigma<- optconstant * zigma
  } else if (i > 500){ 
    XnbarPrev <- Xnbar
    Xnbar <- (i*Xnbar + MC_chain[i, 1:60])/(i+1)
    zigma <- ((i-1)*zigma + tcrossprod(MC_chain[i, 1:60]) + i*tcrossprod(XnbarPrev) - (i+1)*tcrossprod(Xnbar) + epsilon*diag(rep(1,time)))/i
    zigma<- optconstant * zigma
    #print(zigma)
  }   
}

SsMCMC<- numeric(time)
for(i in 1:time){
  SsMCMC[i] = mean(MC_chain[-(1:10000), i])
}
SsMCMC
s

colnames(MC_chain) <- paste("s", 1:60, sep="")

data_df <- as.data.frame(MC_chain)
TrueValues<- s

#Histograms
par(mfrow=c(3, 3))  
for (i in 1:ncol(data_df)) {
  hist(data_df[-(1:10000), i], main = colnames(data_df)[i], xlab ="", col = "white", border = "black")
  abline(v=TrueValues[i], col="red", lwd=2,lty=1)
}

#Traceplots
par(mfrow=c(3, 3))  
for (i in 1:ncol(data_df)) {
  plot(data_df[, i], type = "l", main = colnames(data_df)[i], xlab ="MCMC iterations", ylab = "", col = "purple")
  abline(h=TrueValues[i], col="red", lwd=2,lty=1)
}

