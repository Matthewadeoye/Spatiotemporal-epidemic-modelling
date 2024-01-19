#MCMC for Gamma_12, Gamma_21, and kappa

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
  return ((n - 1)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
  #return (((n - 1)/2) * (log(y) - log(2 * pi)) - ((y/2) * (t(x) %*% z %*% x)))
}

ku<- seq(0, 50, by=0.05)
log_likelihoodku<- lapply(ku, function(x) logIGMRF1(u, x, R))
plot(ku, log_likelihoodku, type="l")

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

num_iteration<- 5000
MC_chain<- matrix(rep(NA,num_iteration*3), nrow=num_iteration, ncol=3)
MC_chain[1,]<- c(runif(1), runif(1), runif(1))
acceptedG12<- 0
acceptedG21<- 0
ACCEPTEDKAPPA<- 0
likelihoodcurrent=loglikelihood(SimulatedData,r,s,u,G(MC_chain[1,1],MC_chain[1,2]),init.density,e.it)
likelihhoodcurrentk<- logIGMRF1(u, MC_chain[1,3], R)

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
  
  proposedkappa<- abs(rnorm(1,mean=MC_chain[i-1,3], sd= 4.8))
  #proposedkappa<- rgamma(1, shape = 1 + 29/2, rate = 1 + t(u) %*% R %*% u /2)
  
  #proposalcurrentkappa<- dgamma(MC_chain[i-1,3], shape = 1 + 29/2, rate = 1 + t(u) %*% R %*% u /2, log=TRUE)
  #proposalproposedkappa<- dgamma(proposedG21, shape = 1 + 29/2, rate = 1 + t(u) %*% R %*% u /2, log=TRUE)

  likelihoodkappa<- logIGMRF1(u, proposedkappa, R)
  likelihoodcurrentk<- logIGMRF1(u, MC_chain[i-1,3], R)
  
  priorcurrentkappa<- dgamma(MC_chain[i-1,3], shape = 1, rate = 0.01, log=TRUE)
  priorproposedkappa<- dgamma(proposedkappa, shape = 1, rate = 0.01, log=TRUE)
  
  #proposalproposedkappa - 
  #+ proposalcurrentkappa
  mh.ratiokappa<- exp(likelihoodkappa + priorproposedkappa - likelihoodcurrentk - priorcurrentkappa)
  print(mh.ratiokappa)
  if(!is.na(mh.ratiokappa) && runif(1) < mh.ratiokappa){
    MC_chain[i,3]<- proposedkappa
    #likelihoodcurrentk=likelihoodkappa
    ACCEPTEDKAPPA<- ACCEPTEDKAPPA + 1
  }
  else{
    MC_chain[i,3]<- MC_chain[i-1,3]
  }
}

#Plots
par(mfrow=c(3,3)) # Create a 2x3 plotting area
hist(MC_chain[500:num_iteration,1], freq=F, xlab = expression(gamma[italic("12")]), breaks=35, main =bquote("Marginal posterior for " * bold(gamma[italic("12")])), col = "white", ylab="Density")
abline(v=0.2, col="red", lwd=2,lty=1)
plot(MC_chain[,1], type = "l", ylab=expression(gamma[italic("12")]), main =bquote("Trace plot for " * bold(gamma[italic("12")])), xlab="MCMC iteration")
plot(MC_chain[, 1], MC_chain[, 2], pch = 16, col = "#00000005", main = "Joint posterior", xlab = expression(gamma[italic("12")]), ylab = expression(gamma[italic("21")]))
hist(MC_chain[500:num_iteration,2], freq=F, xlab = expression(gamma[italic("21")]), breaks=25, main =bquote("Marginal posterior for " * bold(gamma[italic("21")])), col = "white",ylab="Density")
abline(v=0.4, col="red", lwd=2,lty=1)
plot(MC_chain[,2], type = "l", ylab=expression(gamma[italic("21")]), main =bquote("Trace plot for " * bold(gamma[italic("21")])), xlab="MCMC iteration")
plot(MC_chain[,2],MC_chain[,1],  pch = 16, col = "#00000005", main="Joint posterior", xlab= expression(gamma[italic("21")]), ylab=expression(gamma[italic("12")]))
hist(MC_chain[500:num_iteration,3], freq=F, xlab = expression(kappa), breaks=35, main =bquote("Marginal posterior for " * bold(kappa)), col = "white", ylab="Density")
abline(v=kappa_u, col="red", lwd=2,lty=1)
plot(MC_chain[,3], type = "l", ylab=expression(kappa), main =bquote("Trace plot for " * bold(kappa)), xlab="MCMC iteration")

mean(MC_chain[,1])
mean(MC_chain[,2])
mean(MC_chain[,3])


#Frequentist optimization
objective_function <- function(params) {
  
  gammaonetwo<- params[1]
  gammatwoone<- params[2]
  KAPPAU<- params[3]
  
  if (gammaonetwo>1 || gammaonetwo<0 || gammatwoone>1 || gammatwoone<0 ) return(Inf)
  
  m<- matrix(c(1-gammaonetwo,gammaonetwo,gammatwoone,1-gammatwoone),2,2,byrow=T)
  
  llhd<- -(loglikelihood(SimulatedData,r,s,u,Gamma=m,init.density,e.it) + logIGMRF1(u, KAPPAU, R))
  
  return(llhd)
}

initial_params <- c(0.7,0.1, 3)

opt_result <- optim(par = initial_params, fn = objective_function , method = "Nelder-Mead" ,hessian=TRUE)
#print(opt_result)
estimated_params <- opt_result$par
print(estimated_params)



###IGMRF alternative forms
logIGMRF2<- function(x, y) {
  n = nrow(y)
  EV = eigen(Q)$values
  EV = EV[EV>0.0000001]
  return (-(((n - 1)/2) * log(2 * pi)) + (0.5 * log(prod(EV))) - (0.5 * (t(x) %*% y %*% x)))
}

logIGMRF3<- function(u, kapp, node1, node2){
  return((-0.5 * kapp) * (u[node1] - u[node2]) %*% (u[node1] - u[node2]))
}

library(geostan)
stanlist<- prep_icar_data(Westmidlands_adjmat)
N_edges<- stanlist$n_edges
node1<- stanlist$node1
node2<- stanlist$node2

