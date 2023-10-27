# 94 departments  ==> i=1,...,94.
#Monthly cases over 13 years 1985-1997  ==> t=1,...,156.

SimulateMeningococcalData <- function(ndept, time, e_it=matrix(1000,ndept,time)){
  
  #Trend component (r_t)
  
    r <- numeric(time)
    
    kappa_r<- rgamma(1,1,1)
    sigma_r <- sqrt(1/kappa_r)
    
    # initial values
    r[1] <- 0
    r[2] <- 1
    
    # Simulating the trend component
    for (t in 3:time) {
      epsilon <- rnorm(1, mean = 0, sd = sigma_r)  
      r[t] <- (2* (r[t - 1])) - r[t - 2] + epsilon
    }
  
  #Seasonal component (s_t)
  
    s <- numeric(time)
    
    kappa_s<- rgamma(1,1,1)
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

  y_it= matrix(rpois(ndept*time,e_it*exp(logLambda_it)),ndept,time)

return(y_it)

}

SimulatedData<- SimulateMeningococcalData(94,156)

#5th department
plot(SimulatedData[5,],type="l", xlab="Month", ylab="Case")
