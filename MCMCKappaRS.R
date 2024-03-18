#MCMC for Gamma_12, Gamma_21, kappa_r and kappa_s

source("Simulation.R")
source("loglikelihood.R")

SimulatedData<- SimulationResults[[1]]
init.density<- c(0.4, 0.6)
e.it<- 1000

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

kr<- seq(0, 50000, by=0.05)
log_likelihoodkr<- lapply(kr, function(x) randomwalk2(r, x))
plot(kr, log_likelihoodkr, type="l")

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

ks<- seq(0, 5000, by=0.5)
log_likelihoodks<- lapply(ks, function(x) seasonalComp(s, x))
plot(ks, log_likelihoodks, type="l")

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

num_iteration<- 3000
MC_chain<- matrix(NA, nrow=num_iteration, ncol=4)
MC_chain[1,]<- c(runif(1), runif(1), 0, 0)
acceptedG12<- 0
acceptedG21<- 0
ACCEPTEDKAPPA_R<- 0
ACCEPTEDKAPPA_S<- 0
likelihoodcurrent<- loglikelihood(SimulatedData,r,s,u,G(MC_chain[1,1],MC_chain[1,2]),init.density,e.it)

for(i in 2:num_iteration){
  proposedG12<- abs(rnorm(1,mean=MC_chain[i-1,1], sd=0.02))
  if(proposedG12>1) proposedG12=2-proposedG12
  
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
  
  proposedkappaR<- abs(rnorm(1,mean=MC_chain[i-1,3], sd= 13000))
  
  likelihoodproposedkappaR<- randomwalk2(r, proposedkappaR)
  likelihoodcurrentkappaR<- randomwalk2(r, MC_chain[i-1,3])
  
  priorcurrentkappaR<- dgamma(MC_chain[i-1,3], shape = 1, rate = 0.0001, log=TRUE)
  priorproposedkappaR<- dgamma(proposedkappaR, shape = 1, rate = 0.0001, log=TRUE)
  
  mh.ratiokappaR<- exp(likelihoodproposedkappaR + priorproposedkappaR - likelihoodcurrentkappaR - priorcurrentkappaR)
  #print(mh.ratiokappaR)
  if(!is.na(mh.ratiokappaR) && runif(1) < mh.ratiokappaR){
    MC_chain[i,3]<- proposedkappaR
    ACCEPTEDKAPPA_R<- ACCEPTEDKAPPA_R + 1
  }
  else{
    MC_chain[i,3]<- MC_chain[i-1,3]
  }
  
  proposedkappaS<- abs(rnorm(1,mean=MC_chain[i-1,4], sd= 1200))
  
  likelihoodproposedkappaS<- seasonalComp(s, proposedkappaS)
  likelihoodcurrentkappaS<- seasonalComp(s, MC_chain[i-1,4])
  
  priorcurrentkappaS<- dgamma(MC_chain[i-1,4], shape = 1, rate = 0.0001, log=TRUE)
  priorproposedkappaS<- dgamma(proposedkappaS, shape = 1, rate = 0.0001, log=TRUE)
  
  mh.ratiokappaS<- exp(likelihoodproposedkappaS + priorproposedkappaS - likelihoodcurrentkappaS - priorcurrentkappaS)
  #print(mh.ratiokappaS)
  if(!is.na(mh.ratiokappaS) && runif(1) < mh.ratiokappaS){
    MC_chain[i,4]<- proposedkappaS
    ACCEPTEDKAPPA_S<- ACCEPTEDKAPPA_S + 1
  }
  else{
    MC_chain[i,4]<- MC_chain[i-1,4]
  }
}

#Plots
par(mfrow=c(3,3)) # Create a 2x3 plotting area
hist(MC_chain[500:num_iteration,1], freq=F, xlab = expression(gamma[italic("12")]), breaks=35, main =bquote("Marginal posterior for " * bold(gamma[italic("12")])), col = "white", ylab="Density")
abline(v=0.2, col="red", lwd=2,lty=1)
plot(MC_chain[,1], type = "l", ylab=expression(gamma[italic("12")]), main =bquote("Trace plot for " * bold(gamma[italic("12")])), xlab="MCMC iteration")
plot(MC_chain[, 1], MC_chain[, 2], pch = 16, col = "#00000005", main = "Joint posterior FOR Gammas", xlab = expression(gamma[italic("12")]), ylab = expression(gamma[italic("21")]))
hist(MC_chain[500:num_iteration,2], freq=F, xlab = expression(gamma[italic("21")]), breaks=25, main =bquote("Marginal posterior for " * bold(gamma[italic("21")])), col = "white",ylab="Density")
abline(v=0.4, col="red", lwd=2,lty=1)
plot(MC_chain[,2], type = "l", ylab=expression(gamma[italic("21")]), main =bquote("Trace plot for " * bold(gamma[italic("21")])), xlab="MCMC iteration")
#plot(MC_chain[,2],MC_chain[,1],  pch = 16, col = "#00000005", main="Joint posterior", xlab= expression(gamma[italic("21")]), ylab=expression(gamma[italic("12")]))
hist(MC_chain[1000:num_iteration,3], freq=F, xlab = expression(kappa), breaks=35, main =bquote("Marginal posterior for " * bold(kappa)), col = "white", ylab="Density")
abline(v=kappa_r, col="red", lwd=2,lty=1)
plot(MC_chain[,4], type = "l", ylab=expression(kappa), main =bquote("Trace plot for " * bold(kappa)), xlab="MCMC iteration")
hist(MC_chain[1000:num_iteration,4], freq=F, xlab = expression(kappa), breaks=35, main =bquote("Marginal posterior for " * bold(kappa)), col = "white", ylab="Density")
abline(v=kappa_s, col="red", lwd=2,lty=1)
plot(MC_chain[,4], type = "l", ylab=expression(kappa), main =bquote("Trace plot for " * bold(kappa)), xlab="MCMC iteration")

mean(MC_chain[,1])
mean(MC_chain[,2])
mean(MC_chain[,3])
mean(MC_chain[,4])

#Frequentist optimization
objective_function <- function(params) {
  
  gammaonetwo<- params[1]
  gammatwoone<- params[2]
  KAPPAR<- params[3]
  KAPPAS<- params[4]
  
  if (gammaonetwo>1 || gammaonetwo<0 || gammatwoone>1 || gammatwoone<0 ) return(Inf)
  
  m<- matrix(c(1-gammaonetwo,gammaonetwo,gammatwoone,1-gammatwoone),2,2,byrow=T)
  
  llhd<- -(loglikelihood(SimulatedData,r,s,u,Gamma=m,init.density,e.it) + randomwalk2(r, KAPPAR) + seasonalComp(s, KAPPAS))
  
  return(llhd)
}

initial_params <- c(0.7,0.1, 3, 2)

opt_result <- optim(par = initial_params, fn = objective_function , method = "Nelder-Mead" ,hessian=TRUE)
#print(opt_result)
estimated_params <- opt_result$par
print(estimated_params)
