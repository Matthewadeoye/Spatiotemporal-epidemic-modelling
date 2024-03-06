set.seed(122)
source("Simulation.R")
init.density<- c(0.4, 0.6)
stationary_distribution<- c(0.6666667, 0.3333333)
e.it<- 1000

#Function for transition probability matrix
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

#################################################################################################
#FORWARD FILTER
#################################################################################################

forwardfilter <- function(y,r,s,u,Gamma,init.density,e.it) {
  ndept <- nrow(y)  
  nstate <- ncol(Gamma)  
  time <- ncol(y)  
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])
  init.density<- log(init.density)
  
  overallforwardprob<- matrix(0, nrow = ndept, ncol = time)
  
  for (i in 1:ndept) {
    Forwardprob<- matrix(NA, nrow = time, ncol = nstate)
    alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
    alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + 1), log = TRUE)
    
    Forwardprob[1, ] <- c(alpha.1, alpha.2) 
    
    for (t in 2:time) {
      Alphas<- c(alpha.1, alpha.2)
      alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE)))
      alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + 1), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + 1), log = TRUE)))
      Forwardprob[t, ] <- c(alpha.1, alpha.2)
    }
    
    overallforwardprob[i,] <- Forwardprob[, 2]
  }
  
  return(overallforwardprob)
}

#################################################################################################
#BACKWARD SWEEP
#################################################################################################

backwardsweep <- function(y,r,s,u,Gamma,init.density,e.it) {
  ndept <- nrow(y)  
  nstate <- ncol(Gamma)  
  time <- ncol(y)  
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])
  init.density<- log(init.density)
  
overallbackwardprob<- matrix(0, nrow = ndept, ncol = time)
  
  for (i in 1:ndept) {
    Backwardprob<- matrix(NA, nrow = time, ncol = nstate)
    beta.1 <- 0
    beta.2 <- 0
    Backwardprob[time, ] <- c(beta.1, beta.2)
    
    for (t in (time-1):1) {
      Betas<- c(beta.1, beta.2)
      beta.1 <- logSumExp(c(Betas[1] + gamma_11 + dpois(y[i, t+1], lambda = e.it * exp(r[t+1] + s[t+1] + u[i]), log = TRUE), Betas[2] + gamma_21 + dpois(y[i, t+1], lambda = e.it * exp(r[t+1] + s[t+1] + u[i]), log = TRUE)))
      beta.2 <- logSumExp(c(Betas[1] + gamma_12 + dpois(y[i, t+1], lambda = e.it * exp(r[t+1] + s[t+1] + u[i] + 1), log = TRUE), Betas[2] + gamma_22 + dpois(y[i, t+1], lambda = e.it * exp(r[t+1] + s[t+1] + u[i] + 1), log = TRUE)))
      Backwardprob[t, ] <- c(beta.1, beta.2)
    }
    overallbackwardprob[i, ] <- Backwardprob[, 2]
  }
return(overallbackwardprob)
}
