#MCMC for Gamma_12, Gamma_21, and U's
set.seed(122)
source("Simulation.R")
source("loglikelihood.R")

init.density<- c(0.4, 0.6)
stationary_distribution<- c(0.6666667, 0.3333333)
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

num_iteration<- 10000
MC_chain<- matrix(rep(NA,num_iteration*33), nrow=num_iteration, ncol=33)
MC_chain[1,]<- c(runif(1), runif(1), runif(1), rep(0, ndept))
accepted<- 0
likelihoodcurrent<- loglikelihood(SimulatedData,r,s,MC_chain[1,4:33],G(MC_chain[1,1],MC_chain[1,2]),init.density,e.it)
zigma<- diag(rep(0.0001, ndept-1), nrow=ndept-1, ncol=ndept-1)
optconstant<- 2.38^2/(ndept-1)

for(i in 2:num_iteration){
  proposedG12<- abs(rnorm(1,mean=MC_chain[i-1,1], sd=0.02))
  if(proposedG12>1) proposedG12=2-proposedG12
  
  proposedG21<- abs(rnorm(1,mean=MC_chain[i-1,2], sd=0.02))
  if(proposedG21>1) proposedG21=2-proposedG21
  
  proposedkappa<- abs(rnorm(1,mean=MC_chain[i-1,3], sd= 4.8))
  
  proposedspatcomps<- rmvnorm(1, mean=MC_chain[i-1, 4:32], sigma = zigma) 
  proposedspatcomps<- c(proposedspatcomps, -sum(proposedspatcomps))
  
  
  priorcurrent12<- dbeta(MC_chain[i-1,1], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed12<- dbeta(proposedG12, shape1 = 1, shape2 = 1, log=TRUE) 
  
  priorcurrent21<- dbeta(MC_chain[i-1,2], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed21<- dbeta(proposedG21, shape1 = 1, shape2 = 1, log=TRUE)
  
  priorcurrentkappa<- dgamma(MC_chain[i-1,3], shape = 1, rate = 0.01, log=TRUE)
  priorproposedkappa<- dgamma(proposedkappa, shape = 1, rate = 0.01, log=TRUE)
  
  priorcurrentcomps<- logIGMRF1(MC_chain[i-1, 4:33], MC_chain[i-1, 3], R)
  priorproposedcomps<- logIGMRF1(proposedspatcomps, proposedkappa, R)
  
  
  likelihoodproposed<- loglikelihood(SimulatedData,r,s,proposedspatcomps,G(proposedG12,proposedG21),init.density,e.it)

  proposalproposedcomps<- sum(dmvnorm(proposedspatcomps[-30], mean = MC_chain[i-1, 4:32], sigma = zigma, log = TRUE))
  proposalcurrentcomps<- sum(dmvnorm(MC_chain[i-1, 4:32], mean = proposedspatcomps[-30], sigma = zigma, log = TRUE))
  
  mh.ratio<- exp(likelihoodproposed + priorproposed12  + priorproposed21 + priorproposedkappa + priorproposedcomps + proposalcurrentcomps
                  - likelihoodcurrent - priorcurrent12 - priorcurrent21 - priorcurrentkappa - priorcurrentcomps - proposalproposedcomps)
  
  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    accepted<- accepted + 1
    MC_chain[i,1]<- proposedG12
    MC_chain[i,2]<- proposedG21
    MC_chain[i,3]<- proposedkappa
    MC_chain[i,4:33]<- proposedspatcomps
    
    likelihoodcurrent=likelihoodproposed
  }
  else{
    MC_chain[i,1]<- MC_chain[i-1,1]
    MC_chain[i,2]<- MC_chain[i-1,2]
    MC_chain[i,3]<- MC_chain[i-1,3]
    MC_chain[i,4:33]<- MC_chain[i-1,4:33]
    
    accepted<- accepted + 0
  }
  
  #Adapting zigma
  if(i==1500){
    epsilon<- 0.009
    Xn<- MC_chain[1:i, 4:32]
    Xnbar <- colMeans(Xn) 
    zigma <- cov(Xn) + epsilon*diag(rep(1, ndept-1))
    zigma<- optconstant * zigma
  } else if (i > 1500){ 
    XnbarPrev <- Xnbar
    Xnbar <- (i*Xnbar + MC_chain[i, 4:32])/(i+1)
    zigma <- ((i-1)*zigma + tcrossprod(MC_chain[i, 4:32]) + i*tcrossprod(XnbarPrev) - (i+1)*tcrossprod(Xnbar) + epsilon*diag(rep(1,ndept-1)))/i
    zigma<- optconstant * zigma
    print(zigma)
  }   
} 

#Results
UsMCMC<- numeric(ndept)
for(i in 1:ndept){
  UsMCMC[i] = mean(MC_chain[,3+i])
}
sum(UsMCMC)
UsMCMC
u

colnames(MC_chain) <- paste(c("G12", "G21", "kappa_u", paste("u", 1:30, sep="")))

data_df <- as.data.frame(MC_chain)
TrueValues<- c(0.2, 0.4, kappa_u, u)

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

