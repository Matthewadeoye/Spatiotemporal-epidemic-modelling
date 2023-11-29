# 94 departments  ==> i=1,...,94.
#Monthly cases over 13 years 1985-1997  ==> t=1,...,156.

library(ar.matrix)

SimulateMeningococcalData <- function(adj.matrix, time, e_it=matrix(1000,length(diag(adj.matrix)),time)){
  
  ndept<- length(diag(adj.matrix))
  
  #Trend component (r_t)
    #set.seed(8)
    r <- numeric(time)
    
    kappa_r<- rgamma(1,10,0.001)
    sigma_r <- sqrt(1/kappa_r)
    
    # initial values
    r[1] <- runif(1)
    r[2] <- runif(1)
    
    # Simulating the trend component
    for (t in 3:time) {
      epsilon <- rnorm(1, mean = 0, sd = sigma_r)  
      r[t] <- 2*(r[t - 1]) - r[t - 2] + epsilon
    }
  
  #Seasonal component (s_t)
    #set.seed(8)
    s <- numeric(time)
    
    kappa_s<- rgamma(1,10,0.009)
    sigma_s <- sqrt(1/kappa_s)
    
    initial_WN <- rnorm(time, mean = 0, sd = sigma_s)
    
    # Simulating the seasonal component
    for (t in 1:time) {
      if (t <= 11) {
        s[t] <- initial_WN[t]
      } else {
        s[t] <- rnorm(1, mean = -sum(s[(t-1):(t-11)]), sd = sigma_s)
      }
    }
    
  
  #Markov chain (x_it)
  x_it<- function(num_steps){
    
    transition_matrix <- matrix(c(0.8, 0.2, 0.4, 0.6), nrow = 2, byrow = TRUE)
    
    num_steps <- num_steps
    initial_state <- 0  
    
    states <- numeric(num_steps)
    
    # Simulating the Markov chain
    states[1] <- initial_state
    for (i in 2:num_steps) {
      current_state <- states[i - 1]
      next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
      states[i] <- next_state
    }
    return(states)
  }
  
  
  # Using the trend and seasonal components and Markov chain to simulate logLambda_it
  logLambda_it<- matrix(NA, ndept, time)
  for (i in 1:ndept){
    logLambda_it[i,] <- r + s + x_it(time)
  }
  
  #Spatial component (u_i)
  #set.seed(8)
  Q<- adj.matrix
  Q<- -1*Q
  kappa_u<- rgamma(1,10,0.5)
  sigma<- 1/sqrt(kappa_u)
  
  diag(Q)<- -rowSums(Q, na.rm = T)
  Q<- (1/sigma^2)*Q
  
  u<- c(sim.GMRF(1,Q))
  
  logLambda_it = logLambda_it + u

  y_it= matrix(rpois(ndept*time,e_it*exp(logLambda_it)),ndept,time)

return(y_it)

}

#Usage

#ADJACENCY MATRIX FRENCH DEPARTMENTS
library(sf)
library(mapview)
library(spdep)
Westmidlands.shp <- st_read("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/LAD_MAY_2023_UK_BUC_V2.shp")
mapview(Westmidlands.shp)
Westmidlands_adjmat <- nb2mat(poly2nb(Westmidlands.shp), style = "B", zero.policy = TRUE)
Westmidlands_adjmat[1:10, 1:10]
adj.matrix<- Westmidlands_adjmat

set.seed(8)
SimulatedData<- SimulateMeningococcalData(adj.matrix=Westmidlands_adjmat, time=24)

#5th department
plot(SimulatedData[5,],type="l", xlab="Month", ylab="Case")


#Likelihood Computation via forward filtering
y <- SimulatedData
Gamma <- matrix(c(0.8, 0.2, 0.4, 0.6), nrow = 2, byrow = TRUE)
init.density <- c(0.5, 0.5)
ndept <- nrow(y)  
time <- ncol(y)  
nstate <- ncol(Gamma)  
e.it<- 1000

likelihood <- function(y,r,s,u,Gamma,init.density,e.it) {
  forward.probs <- matrix(NA, ndept, nstate)
  for (i in 1:ndept) {
    alpha <- init.density * c(dpois(y[i, 1], lambda = e.it*(exp(r[1] + s[1] + u[i]))),
                              dpois(y[i, 1], lambda = e.it*(exp(r[1] + s[1] + u[i] + 1))))
    logAlph<- log(sum(alpha))
    
    for (t in 2:time) {
      alpha <- alpha %*% Gamma * c(dpois(y[i, t], lambda = e.it*(exp(r[t] + s[t] + u[i]))),
                                   dpois(y[i, t], lambda = e.it*(exp(r[t] + s[t] + u[i] + 1))))
      logAlph<- logAlph + log(sum(alpha))
    }
    forward.probs[i,] <- alpha
  }
  lklhd<- forward.probs[,1] + forward.probs[,2]
  full.likelihood<- prod(lklhd)
  full.log.likelihood<- logAlph
  return(list(forward.probs=forward.probs, full.likelihood=full.likelihood, full.log.likelihood=full.log.likelihood))
}

likelihood(y,r,s,u,Gamma,init.density,e.it)


#########################################################################################################

#logLikelihood Computation via forward filtering
loglikelihood <- function(y,r,s,u,Gamma,init.density,e.it) {
  ndept <- nrow(y)  
  nstate <- ncol(Gamma)  
  time <- ncol(y)  
  init.density<- log(init.density)
  
  log.forward.probs <- matrix(NA, ndept, nstate)
  
  for (i in 1:ndept) {
    
    alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
    alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + 1), log = TRUE)
    
    for (t in 2:time) {
      Alphas<- c(alpha.1, alpha.2)
      if (Alphas[2] + log(Gamma[2,1]) > Alphas[1] + log(Gamma[1,1])) 
        alpha.1 <- Alphas[2] + log(Gamma[2,1]) + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE) + log(1 + exp(Alphas[1] + log(Gamma[1,1]) - Alphas[2] - log(Gamma[2,1])))
        else
          alpha.1 <- Alphas[1] + log(Gamma[1,1]) + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE) + log(1 + exp(Alphas[2] + log(Gamma[2,1]) - Alphas[1] - log(Gamma[1,1])))
    
            if(Alphas[2] + log(Gamma[2,2]) > Alphas[1] + log(Gamma[1,2])) 
          alpha.2 <- Alphas[2] + log(Gamma[2,2]) + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + 1), log = TRUE) + log(1 + exp(Alphas[1] + log(Gamma[1,2]) - Alphas[2] - log(Gamma[2,2])))
          else
           alpha.2 <- Alphas[1] + log(Gamma[1,2]) + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + 1), log = TRUE) + log(1 + exp(Alphas[2] + log(Gamma[2,2]) - Alphas[1] - log(Gamma[1,2])))
    }
    
    log.forward.probs[i,] <- c(alpha.1, alpha.2)
  }
  
  full.log.likelihood<- 0
  for (i in 1:ndept) {
  if (log.forward.probs[i,1] > log.forward.probs[i,2]) {
    M <- log.forward.probs[i,1]
  } else {
    M <- log.forward.probs[i,2]
  }
  full.log.likelihood<- full.log.likelihood + M + log(1+ exp(min(log.forward.probs[i,1],log.forward.probs[i,2]) - M))
  }
  
  return(list(log.forward.probs=log.forward.probs, full.log.likelihood=full.log.likelihood))
}

loglikelihood(y,r,s,u,Gamma,init.density,e.it)

