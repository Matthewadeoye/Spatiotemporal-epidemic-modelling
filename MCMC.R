#MCMC for Gamma_12 and Gamma_21

source("Simulation.R")
source("loglikelihood.R")

init.density<- c(0.4, 0.6)
stationary_distribution<- c(0.6666667, 0.3333333)
e.it<- 1000

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

num_iteration<- 5000
MC_chain<- matrix(rep(NA,num_iteration*2), nrow=num_iteration, ncol=2)
MC_chain[1,]<- c(runif(1), runif(1))
acceptedG12<- 0
acceptedG21<- 0
likelihoodcurrent=loglikelihood(SimulatedData,r,s,u,G(MC_chain[1,1],MC_chain[1,2]),init.density,e.it)

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
}

#Plots
par(mfrow=c(2,3)) # Create a 2x3 plotting area
hist(MC_chain[,1], freq=F, xlab = expression(gamma[italic("12")]), breaks=35, main =bquote("Marginal posterior for " * bold(gamma[italic("12")])), col = "white", ylab="Density")
abline(v=0.2, col="red", lwd=2,lty=1)
plot(MC_chain[,1], type = "l", ylab=expression(gamma[italic("12")]), main =bquote("Trace plot for " * bold(gamma[italic("12")])), xlab="MCMC iteration")
plot(MC_chain[, 1], MC_chain[, 2], pch = 16, col = "#00000005", main = "Joint posterior", xlab = expression(gamma[italic("12")]), ylab = expression(gamma[italic("21")]))
hist(MC_chain[,2], freq=F, xlab = expression(gamma[italic("21")]), breaks=25, main =bquote("Marginal posterior for " * bold(gamma[italic("21")])), col = "white",ylab="Density")
abline(v=0.4, col="red", lwd=2,lty=1)
plot(MC_chain[,2], type = "l", ylab=expression(gamma[italic("21")]), main =bquote("Trace plot for " * bold(gamma[italic("21")])), xlab="MCMC iteration")
plot(MC_chain[,2],MC_chain[,1],  pch = 16, col = "#00000005", main="Joint posterior", xlab= expression(gamma[italic("21")]), ylab=expression(gamma[italic("12")]))

mean(MC_chain[,1])
mean(MC_chain[,2])
#acceptedG12/num_iteration
#acceptedG21/num_iteration



#Frequentist optimization
objective_function <- function(params) {
  
  gammaonetwo<- params[1]
  gammatwoone<- params[2]

  if (gammaonetwo>1 || gammaonetwo<0 || gammatwoone>1 || gammatwoone<0 ) return(Inf)
  
  m<- matrix(c(1-gammaonetwo,gammaonetwo,gammatwoone,1-gammatwoone),2,2,byrow=T)
  
  llhd<- -loglikelihood(SimulatedData,r,s,u,Gamma=m,init.density,e.it)
  
  return(llhd)
}

initial_params <- c(0.7,0.1)

opt_result <- optim(par = initial_params, fn = objective_function , method = "CG" ,hessian=TRUE)
#print(opt_result)
estimated_params <- opt_result$par
print(estimated_params)



#stationary distribution
TPM <- G(0.2, 0.4)

eigen_decomp <- eigen(t(TPM))
eigenvalues <- eigen_decomp$values
eigenvectors <- eigen_decomp$vectors

index <- which(abs(eigenvalues - 1) < 1e-8)

stationary_distribution <- eigenvectors[, index]

stationary_distribution <- stationary_distribution / sum(stationary_distribution)

stationary_distribution
