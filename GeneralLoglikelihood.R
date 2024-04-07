#logLikelihood via forward filtering

GeneralLoglikelihood <- function(y,r,s,u,Gamma,init.density,e.it,B,model,adjacencyMatrix) {
  ndept <- nrow(y)  
  nstate <- ncol(Gamma)  
  time <- ncol(y)  
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])
  init.density<- log(init.density)
  
  log.forward.probs <- matrix(NA, ndept, nstate)
  
  if(model == 0){
  
  for (i in 1:ndept) {
    alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
    alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + 1), log = TRUE)
    
    for (t in 2:time) {
      Alphas<- c(alpha.1, alpha.2)
      alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE)))
      alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + 1), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + 1), log = TRUE)))
    }
    
    log.forward.probs[i,] <- c(alpha.1, alpha.2)
  }
  
  full.log.likelihood<- rowLogSumExps(log.forward.probs)
  
  return(sum(full.log.likelihood))
 }
 else
  if(model == 1){
    #Model1
    z.it<- matrix(0, ndept, time)
    for(i in 1:ndept){
      for(t in 2:time){
        if(y[i, t-1] > 0)
          z.it[i, t]<- 1
      }
    }
    
    for (i in 1:ndept) {
      alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
      alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + B[1] * z.it[i,1]), log = TRUE)
      
      for (t in 2:time) {
        Alphas<- c(alpha.1, alpha.2)
        alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE)))
        alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * z.it[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * z.it[i,t]), log = TRUE)))
      }
      
      log.forward.probs[i,] <- c(alpha.1, alpha.2)
    }
    
    full.log.likelihood<- rowLogSumExps(log.forward.probs)
    
    return(sum(full.log.likelihood))
  }
  else
    if(model == 2){
      #Model2
      z.it<- matrix(0, ndept, time)
      
      for(i in 1:ndept){
        
        indexes<- c()
        for(a in 1:ndept){
          if(i != a){
            if(adjacencyMatrix[i,a] > 0)
              indexes<- c(indexes, a)
          }
        }
        
        for(t in 2:time){
          for(b in 1:length(indexes)){
          if(y[i, t-1] > 0 | y[indexes[b], t-1] > 0)
            z.it[i, t]<- 1
          }
        }
      }
        
      for (i in 1:ndept) {
        alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + B[1] * z.it[i,1]), log = TRUE)
        
        for (t in 2:time) {
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * z.it[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * z.it[i,t]), log = TRUE)))
        }
        
        log.forward.probs[i,] <- c(alpha.1, alpha.2)
      }
      
      full.log.likelihood<- rowLogSumExps(log.forward.probs)
      
      return(sum(full.log.likelihood))
    }
  else
    if(model == 3){
      #Model3
      z.it<- matrix(0, ndept, time)
      for(i in 1:ndept){
        for(t in 2:time){
          if(y[i, t-1] > 0)
            z.it[i, t]<- 1
        }
      }
      
      z.it2<- matrix(0, ndept, time)
      for(i in 1:ndept){
        
        indexes<- c()
        for(a in 1:ndept){
          if(i != a){
            if(adjacencyMatrix[i,a] > 0)
              indexes<- c(indexes, a)
          }
        }
        
        for(t in 2:time){
          for(b in 1:length(indexes)){
            if(y[indexes[b], t-1] > 0)
              z.it2[i, t]<- 1
          }
        }
      }
      
      for (i in 1:ndept) {
        alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + B[1] * z.it[i,1] + B[2] * z.it2[i,1]), log = TRUE)
        
        for (t in 2:time) {
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * z.it[i,t] + B[2] * z.it2[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * z.it[i,t] + B[2] * z.it2[i,t]), log = TRUE)))
        }
        
        log.forward.probs[i,] <- c(alpha.1, alpha.2)
      }
      
      full.log.likelihood<- rowLogSumExps(log.forward.probs)
      
      return(sum(full.log.likelihood))
    }
  else
    if(model == 4){
      #Model4
      
      for (i in 1:ndept) {
        alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + B[1] * log(0 + 1)), log = TRUE)
        
        for (t in 2:time) {
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + 1)), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + 1)), log = TRUE)))
        }
        
        log.forward.probs[i,] <- c(alpha.1, alpha.2)
      }
      
      full.log.likelihood<- rowLogSumExps(log.forward.probs)
      
      return(sum(full.log.likelihood))
    } 
  else
    if(model == 5){
      #Model5
      z.it<- matrix(0, ndept, time)
      for(i in 1:ndept){
        
        indexes<- c()
        for(a in 1:ndept){
          if(i != a){
            if(adjacencyMatrix[i,a] > 0)
              indexes<- c(indexes, a)
          }
        }
        
        for(t in 2:time){
          for(b in 1:length(indexes)){
              z.it[i, t]<- z.it[i, t] + y[indexes[b], t-1]
          }
        }
      }
      
      for (i in 1:ndept) {
        alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + B[1] * log(0 + 0 + 1)), log = TRUE)
        
        for (t in 2:time) {
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + z.it[i, t] + 1)), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + z.it[i, t] + 1)), log = TRUE)))
        }
        
        log.forward.probs[i,] <- c(alpha.1, alpha.2)
      }
      
      full.log.likelihood<- rowLogSumExps(log.forward.probs)
      
      return(sum(full.log.likelihood))
    } 
  else
    if(model == 6){
      #Model6
      z.it<- matrix(0, ndept, time)
      for(i in 1:ndept){
        
        indexes<- c()
        for(a in 1:ndept){
          if(i != a){
            if(adjacencyMatrix[i,a] > 0)
              indexes<- c(indexes, a)
          }
        }
        
        for(t in 2:time){
          for(b in 1:length(indexes)){
            z.it[i, t]<- z.it[i, t] + y[indexes[b], t-1]
          }
        }
      }
      
      for (i in 1:ndept) {
        alpha.1 <- init.density[1] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init.density[2] + dpois(y[i, 1], lambda = e.it * exp(r[1] + s[1] + u[i] + B[1] * log(0 + 1) + B[2] * z.it[i, 1] + 1), log = TRUE)
        
        for (t in 2:time) {
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + 1) + B[2] * z.it[i, t] + 1), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e.it * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + 1) + B[2] * z.it[i, t] + 1), log = TRUE)))
        }
        
        log.forward.probs[i,] <- c(alpha.1, alpha.2)
      }
      
      full.log.likelihood<- rowLogSumExps(log.forward.probs)
      
      return(sum(full.log.likelihood))
    } 
}


