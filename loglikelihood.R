#logLikelihood via forward filtering

loglikelihood <- function(y,r,s,u,Gamma,init.density,e_it) {
  ndept <- nrow(y)  
  nstate <- ncol(Gamma)  
  time <- ncol(y)  
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])
  init.density<- log(init.density)
  
  log.forward.probs <- matrix(NA, ndept, nstate)
  
  for (i in 1:ndept) {
    alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
    alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + 1), log = TRUE)
    
    for (t in 2:time) {
      Alphas<- c(alpha.1, alpha.2)
      alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE)))
      alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + 1), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + 1), log = TRUE)))
    }
    
    log.forward.probs[i,] <- c(alpha.1, alpha.2)
  }
  
  full.log.likelihood<- rowLogSumExps(log.forward.probs)
  
  return(sum(full.log.likelihood))
}


#logLikelihood0 for model 7

loglikelihood0<- function(y,r,s,u,e_it) {
  ndept <- nrow(y)  
  time <- ncol(y)  
  
  log.likelihood <- matrix(NA, ndept, time)
  
  for (i in 1:ndept) {
    for (t in 2:time) {
      log.likelihood[i, t]<- dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE)
    }
  }
  
  full.log.likelihood<- sum(log.likelihood)
  
  return(full.log.likelihood)
}

#logLikelihood via forward filtering

loglikelihood01 <- function(y,r,s,u,init.density,e_it) {
  ndept <- nrow(y)  
  nstate <- length(init.density)
  time <- ncol(y)  
  gamma_11 <- log(0.9999)
  gamma_12 <- log(0.0001)
  gamma_21 <- log(0.0001)
  gamma_22 <- log(0.9999)
  init.density<- log(init.density)
  
  log.forward.probs <- matrix(NA, ndept, nstate)
  
  for (i in 1:ndept) {
    alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
    alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + 1), log = TRUE)
    
    for (t in 2:time) {
      Alphas<- c(alpha.1, alpha.2)
      alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE)))
      alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + 1), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + 1), log = TRUE)))
    }
    
    log.forward.probs[i,] <- c(alpha.1, alpha.2)
  }
  
  full.log.likelihood<- rowLogSumExps(log.forward.probs)
  
  return(sum(full.log.likelihood))
}