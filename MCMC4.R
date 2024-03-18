#MCMC for Gamma_12, Gamma_21, and U's
set.seed(122)
source("Simulation.R")
source("loglikelihood.R")

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
RW2PrecMat<- (1/kappa_r)*RW2PrecMat

#Prior density for Trend components (r_t)
randomwalk2<- function(componentR, PrecisionR){
  time<- length(componentR)
  Sumres<- 0
  for(i in 3:time){
    res<- (componentR[i-2] - (2 * componentR[i-1]) + componentR[i])^2
    Sumres<- Sumres + res
  }
  return((time - 2)/2 * (log(PrecisionR) - log(2 * pi)) -PrecisionR/2 * Sumres)
}

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

num_iteration<- 4000
MC_chain<- matrix(NA, nrow=num_iteration, ncol=time)
MC_chain[1,]<- r
ACCEPTED<- 0
likelihoodcurrent<- loglikelihood(SimulatedData,MC_chain[1, 1:60],s,u,G(0.2, 0.4),init.density,e.it)
zigma<- RW2PrecMat
diag(zigma)<- diag(zigma) + 0.000001
zigma<- solve(zigma)
optconstant<- 2.38^2/(time)

for(i in 2:num_iteration){
  proposedRcomps<- rmvnorm(1, mean=MC_chain[i-1, 1:60], sigma = zigma) 
  
  likelihoodproposed<- loglikelihood(SimulatedData,proposedRcomps,s,u,G(0.2, 0.4),init.density,e.it)
  
  priorcurrentRcomps<- randomwalk2(MC_chain[i-1, 1:60], kappa_r)
  priorproposedRcomps<- randomwalk2(proposedRcomps, kappa_r)
  
  proposalproposedRcomps<- sum(dmvnorm(proposedRcomps, mean = MC_chain[i-1, 1:60], sigma = zigma, log = TRUE))
  proposalcurrentRcomps<- sum(dmvnorm(MC_chain[i-1, 1:60], mean = proposedRcomps, sigma = zigma, log = TRUE))
  
  
  mh.ratio<- exp(likelihoodproposed + priorproposedRcomps + proposalcurrentRcomps - likelihoodcurrent - priorcurrentRcomps - proposalproposedRcomps)
  
  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    ACCEPTED<- ACCEPTED + 1
    MC_chain[i,1:60]<- proposedRcomps
    likelihoodcurrent<- likelihoodproposed
  }
  else{
    MC_chain[i, 1:60]<- MC_chain[i-1, 1:60]
    ACCEPTED<- ACCEPTED + 0
  }
  
  #Adapting zigma
  if(i==800){
    epsilon<- 0.009
    Xn<- MC_chain[1:i, 1:60]
    Xnbar <- colMeans(Xn) 
    zigma <- cov(Xn) + epsilon*diag(rep(1, time))
    zigma<- optconstant * zigma
  } else if (i > 800){ 
    XnbarPrev <- Xnbar
    Xnbar <- (i*Xnbar + MC_chain[i, 1:60])/(i+1)
    zigma <- ((i-1)*zigma + tcrossprod(MC_chain[i, 1:60]) + i*tcrossprod(XnbarPrev) - (i+1)*tcrossprod(Xnbar) + epsilon*diag(rep(1,time)))/i
    zigma<- optconstant * zigma
    print(zigma)
  }   
}

RsMCMC<- numeric(time)
for(i in 1:time){
  RsMCMC[i] = mean(MC_chain[, i])
}
sum(RsMCMC)
RsMCMC
r

colnames(MC_chain) <- paste("r", 1:60, sep="")

data_df <- as.data.frame(MC_chain)
TrueValues<- r

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

