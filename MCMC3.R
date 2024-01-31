#MCMC for Gamma_12, Gamma_21, and U's
set.seed(122)
source("Simulation.R")
source("loglikelihood.R")

init.density<- c(0.4, 0.6)
e.it<- 1000

#Intrinsic GMRF density
R<- -1*Westmidlands_adjmat
diag(R)<- -rowSums(R, na.rm = T)
qr(R)$rank

#x are the components(u), y is the parameter (kappa_u), z is the structure matrix (R)
logIGMRF1<- function(x, y, z) {
  n = nrow(z)
  sumC = sum(x[1:n-1])
  x = c(x[1:n-1], -sumC)
  return ((n - 1)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
}

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

num_iteration<- 800
MC_chain<- matrix(rep(NA,num_iteration*32), nrow=num_iteration, ncol=32)
MC_chain[1,]<- c(runif(1), runif(1), c(rep(0, ndept)))
acceptedG12<- 0
acceptedG21<- 0
ACCEPTEDKAPPA<- 0
acceptedcomps<- 0
likelihoodcurrent=loglikelihood(SimulatedData,r,s,MC_chain[1,3:32],G(MC_chain[1,1],MC_chain[1,2]),init.density,e.it)
#likelihhoodcurrentk<- logIGMRF1(u, MC_chain[1,3], R)

for(i in 2:num_iteration){
  proposedG12<- abs(rnorm(1,mean=MC_chain[i-1,1], sd=0.02))
  if(proposedG12>1) proposedG12=2-proposedG12
  
  #likelihoodcurrent12<- loglikelihood(SimulatedData,r,s,u,G(MC_chain[i-1,1],MC_chain[i-1,2]),init.density,e.it)
  likelihoodproposed12<- loglikelihood(SimulatedData,r,s,u,G(proposedG12,MC_chain[i-1,2]),init.density,e.it)
  
  priorcurrent12<- dbeta(MC_chain[i-1,1], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed12<- dbeta(proposedG12, shape1 = 1, shape2 = 1, log=TRUE) 
  
  mh.ratioG12<- exp(likelihoodproposed12 + priorproposed12 - likelihoodcurrent - priorproposed12)
  
  if(!is.na(mh.ratioG12) && runif(1) < mh.ratioG12){
    acceptedG12<- acceptedG12 + 1
    MC_chain[i,1]<- proposedG12
    likelihoodcurrent=likelihoodproposed12
  }
  else{
    MC_chain[i,1]<- MC_chain[i-1,1]
    acceptedG12<- acceptedG12 + 0
  }
  
  proposedG21<- abs(rnorm(1,mean=MC_chain[i-1,2], sd=0.02))
  if(proposedG21>1) proposedG21=2-proposedG21
  
  #likelihoodcurrent21<- loglikelihood(SimulatedData,r,s,u,G(MC_chain[i-1,1],MC_chain[i-1,2]),init.density,e.it) 
  likelihoodproposed21<- loglikelihood(SimulatedData,r,s,u,G(MC_chain[i-1,1], proposedG21),init.density,e.it)
  
  priorcurrent21<- dbeta(MC_chain[i-1,2], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed21<- dbeta(proposedG21, shape1 = 1, shape2 = 1, log=TRUE)
  
  
  mh.ratioG21<- exp(likelihoodproposed21 + priorproposed21 - likelihoodcurrent - priorproposed21)
  
  if(!is.na(mh.ratioG21) && runif(1) < mh.ratioG21){
    acceptedG21<- acceptedG21 + 1
    MC_chain[i,2]<- proposedG21
    likelihoodcurrent=likelihoodproposed21
  }
  else{
    MC_chain[i,2]<- MC_chain[i-1,2]
    acceptedG21<- acceptedG21 + 0
  }
  
  #proposedkappa<- abs(rnorm(1,mean=MC_chain[i-1,3], sd= 4.8))
  #proposedkappa<- rgamma(1, shape = 1 + 29/2, rate = 1 + t(u) %*% R %*% u /2)
  
  #proposalcurrentkappa<- dgamma(MC_chain[i-1,3], shape = 1 + 29/2, rate = 1 + t(u) %*% R %*% u /2, log=TRUE)
  #proposalproposedkappa<- dgamma(proposedG21, shape = 1 + 29/2, rate = 1 + t(u) %*% R %*% u /2, log=TRUE)
  
  #likelihoodkappa<- logIGMRF1(u, proposedkappa, R)
  #likelihoodcurrentk<- logIGMRF1(u, MC_chain[i-1,3], R)
  
  #priorcurrentkappa<- dgamma(MC_chain[i-1,3], shape = 1, rate = 0.01, log=TRUE)
  #priorproposedkappa<- dgamma(proposedkappa, shape = 1, rate = 0.01, log=TRUE)
  
  #proposalproposedkappa - 
  #+ proposalcurrentkappa
  #mh.ratiokappa<- exp(likelihoodkappa + priorproposedkappa - likelihoodcurrentk - priorcurrentkappa)
  #print(mh.ratiokappa)
  #if(!is.na(mh.ratiokappa) && runif(1) < mh.ratiokappa){
  #  MC_chain[i,3]<- proposedkappa
    #likelihoodcurrentk=likelihoodkappa
  #  ACCEPTEDKAPPA<- ACCEPTEDKAPPA + 1
  #}
  #else{
  #  MC_chain[i,3]<- MC_chain[i-1,3]
  #}
  
  tunedL<- 3.009
  R<- -1*Westmidlands_adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  diag(R) <- diag(R) + 19.01
  
  zigma<- diag(rep(0.0001, ndept), nrow=ndept, ncol=ndept)
  #zigma<- solve(tunedL*kappa_u*R)
  
  proposedspatcomps<- rmvnorm(1, mean=MC_chain[i-1, 3:32], sigma = zigma) 
  #proposedspatcomps<- c(sim.GMRF(1, tunedL * kappa_u * R))
  
  likelihoodproposedcomps<- loglikelihood(SimulatedData,r,s,proposedspatcomps,G(0.2,0.4),init.density,e.it)
  likelihoodcurrentcomps<- loglikelihood(SimulatedData,r,s,MC_chain[i-1, 3:32],G(0.2,0.4),init.density,e.it)
  
  
  priorcurrentcomps<- logIGMRF1(MC_chain[i-1, 3:32], kappa_u, R)
  priorproposedcomps<- logIGMRF1(proposedspatcomps, kappa_u, R)
  
  proposalproposedcomps<- sum(dmvnorm(proposedspatcomps, mean = MC_chain[i-1, 3:32], sigma = zigma, log = TRUE))
  proposalcurrentcomps<- sum(dmvnorm(MC_chain[i-1, 3:32], mean = proposedspatcomps, sigma = zigma, log = TRUE))
  
  mh.ratiokappa<- exp(likelihoodproposedcomps + priorproposedcomps + proposalcurrentcomps - likelihoodcurrentcomps - priorcurrentcomps - proposalproposedcomps)
  print(mh.ratiokappa)
  
  if(!is.na(mh.ratiokappa) && runif(1) < mh.ratiokappa){
    MC_chain[i,3:32]<- proposedspatcomps
    acceptedcomps<- acceptedcomps + 1
  }
  else{
    MC_chain[i,3:32]<- MC_chain[i-1,3:32]
  }
}

mean(MC_chain[,1])
mean(MC_chain[,2])
mean(MC_chain[,3])

UsMCMC<- numeric(ndept)
for(i in 1:ndept){
  UsMCMC[i] = mean(MC_chain[,2+i])
}
sum(UsMCMC)
UsMCMC
u
plot(MC_chain[,3], type = "l", ylab="chain", main ="Trace plot", xlab="MCMC iteration")

############################################################################################

num_iteration<- 5000
MC_chain<- matrix(NA, nrow=num_iteration, ncol=30)
MC_chain[1,]<- c(rep(0, ndept))
acceptedcomps<- 0

for(i in 2:num_iteration){
tunedL<- 10.009
R<- -1*Westmidlands_adjmat
diag(R)<- -rowSums(R, na.rm = T)
diag(R) <- diag(R) + 16.01

#zigma<- tunedL*kappa_u*R
zigma<- solve(tunedL*kappa_u*R)

proposedspatcomps<- rmvnorm(1, mean=MC_chain[i-1, 1:30], sigma = zigma) 
#proposedspatcomps<- c(sim.GMRF(1, tunedL * kappa_u * R))

likelihoodproposedcomps<- loglikelihood(SimulatedData,r,s,proposedspatcomps,G(0.2,0.4),init.density,e.it)
likelihoodcurrentcomps<- loglikelihood(SimulatedData,r,s,MC_chain[i-1, 1:30],G(0.2,0.4),init.density,e.it)


priorcurrentcomps<- logIGMRF1(MC_chain[i-1, 1:30], kappa_u, R)
priorproposedcomps<- logIGMRF1(proposedspatcomps, kappa_u, R)

proposalproposedcomps<- sum(dmvnorm(proposedspatcomps, mean = MC_chain[i-1, 1:30], sigma = zigma, log = TRUE))
proposalcurrentcomps<- sum(dmvnorm(MC_chain[i-1, 1:30], mean = proposedspatcomps, sigma = zigma, log = TRUE))

mh.ratiokappa<- exp(likelihoodproposedcomps + priorproposedcomps + proposalcurrentcomps - likelihoodcurrentcomps - priorcurrentcomps - proposalproposedcomps)
print(mh.ratiokappa)

if(!is.na(mh.ratiokappa) && runif(1) < mh.ratiokappa){
  MC_chain[i,1:30]<- proposedspatcomps
  acceptedcomps<- acceptedcomps + 1
}
else{
  MC_chain[i,1:30]<- MC_chain[i-1,1:30]
 }
}
plot(MC_chain[,30], type = "l", ylab=expression(kappa), main ="Trace plot", xlab="MCMC iteration")

UsMCMC<- numeric(ndept)
for(i in 1:ndept){
  UsMCMC[i] = mean(MC_chain[,i])
}
sum(UsMCMC)
UsMCMC
u