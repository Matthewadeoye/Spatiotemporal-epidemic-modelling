#logLikelihood via forward filtering

GeneralLoglikelihood <- function(y,r,s,u,Gamma,init_density,e_it,B,model,adjacencyMatrix, z_it, z_it2) {
  ndept <- nrow(y)  
  nstate <- ncol(Gamma)  
  time <- ncol(y)  
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])
  init_density<- log(init_density)
  
  log.forward.probs <- matrix(NA, ndept, nstate)
  
  #Model0
  if(model == 0){
    log.likelihood <- matrix(NA, ndept, time)
    
    for (i in 1:ndept) {
      for (t in 1:time) {
        log.likelihood[i, t]<- dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE)
      }
    }
    
    full.log.likelihood<- sum(log.likelihood)
    
    return(full.log.likelihood)
  }
  
  #Model 1, 2, 4 or 5
 else if(model == 1 || model ==2 || model == 4 || model ==5){
    for (i in 1:ndept) {
      alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
      alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i,1]), log = TRUE)
      
      for (t in 2:time) {
        Alphas<- c(alpha.1, alpha.2)
        alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE)))
        alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i,t]), log = TRUE)))
      }
      
      log.forward.probs[i,] <- c(alpha.1, alpha.2)
    }
    
    full.log.likelihood<- rowLogSumExps(log.forward.probs)
    
    return(sum(full.log.likelihood))
  }
  
  #Model3 or Model6
  else if(model == 3 || model == 6){
      for (i in 1:ndept) {
        alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i,1] + B[2] * z_it2[i,1]), log = TRUE)
        
        for (t in 2:time) {
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i, t] + B[2] * z_it2[i, t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i, t] + B[2] * z_it2[i, t]), log = TRUE)))
        }
        
        log.forward.probs[i,] <- c(alpha.1, alpha.2)
      }
      
      full.log.likelihood<- rowLogSumExps(log.forward.probs)
      
      return(sum(full.log.likelihood))
    }
 
  #Model 7
  else if(model == 7){
      
      for (i in 1:ndept) {
      alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
      alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + 1), log = TRUE)
      
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
}


