#Computing the probabilities of epidemic outbreak

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
  
  AllForwardprobs<- vector("list", ndept)

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
    AllForwardprobs[[i]]<- Forwardprob
  }
  
  return(AllForwardprobs)
}

#################################################################################################
#BACKWARD SWEEP
#################################################################################################

backwardsweep <- function(y,r,s,u,Gamma,e.it) {
  ndept <- nrow(y)  
  nstate <- ncol(Gamma)  
  time <- ncol(y)  
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])

  Allbackwardprob<- vector("list", ndept)
  
  for (i in 1:ndept) {
    Backwardprob<- matrix(NA, nrow = time, ncol = nstate)
    beta.1 <- 0
    beta.2 <- 0
    Backwardprob[time, ] <- c(beta.1, beta.2)
    
    for (t in (time-1):1) {
      Betas<- c(beta.1, beta.2)
      beta.1 <- logSumExp(c(Betas[1] + gamma_11 + dpois(y[i, t+1], lambda = e.it * exp(r[t+1] + s[t+1] + u[i]), log = TRUE), Betas[2] + gamma_12 + dpois(y[i, t+1], lambda = e.it * exp(r[t+1] + s[t+1] + u[i]), log = TRUE)))
      beta.2 <- logSumExp(c(Betas[1] + gamma_21 + dpois(y[i, t+1], lambda = e.it * exp(r[t+1] + s[t+1] + u[i] + 1), log = TRUE), Betas[2] + gamma_22 + dpois(y[i, t+1], lambda = e.it * exp(r[t+1] + s[t+1] + u[i] + 1), log = TRUE)))
      Backwardprob[t, ] <- c(beta.1, beta.2)
    }
    Allbackwardprob[[i]] <- Backwardprob
  }
return(Allbackwardprob)
}


#Local decoding
Decoding <- function(y,r,s,u,Gamma,init.density,e.it) {
  ndept<- nrow(y)
  time <- ncol(y)  
  Allforwardprobs<- forwardfilter(y,r,s,u,Gamma,init.density,e.it)
  Allbackwardprobs<- backwardsweep(y,r,s,u,Gamma,e.it)
  Res<- matrix(NA, ndept, time)  
  for(i in 1:ndept){
    for(j in 1:time){
      P1<- Allforwardprobs[[i]][j,1] + Allbackwardprobs[[i]][j,1]
      P2<- Allforwardprobs[[i]][j,2] + Allbackwardprobs[[i]][j,2]
      Res[i,j]<- exp(P2 - logSumExp(c(P1,P2)))
      #browser()
    }
  }
  return(Res)
}  
  