# 94 departments  ==> i=1,...,94.
# Monthly cases over 13 years 1985-1997  ==> t=1,...,156.

SimulateMeningococcalData <- function(ndept, time){
  time1<- 1:time
  logLambda_it <- c() 
  Month <- rep(time1, ndept)
  Department <- as.factor(rep(1:ndept, each = length(time1)))
  
  #Trend component (r_t)
  r_t<- function(num_steps){
    num_steps <- num_steps
    
    kappa_r<- rgamma(1,1,0.0001)
    sigma_r <- sqrt(1/kappa_r)
    
    r <- numeric(num_steps)
    
    # initial values
    r[1] <- 0
    r[2] <- 1
    
    # Simulating the trend component
    for (t in 3:num_steps) {
      epsilon <- rnorm(1, mean = 0, sd = sigma_r)  
      r[t] <- (2* (r[t - 1])) - r[t - 2] + epsilon
    }
    
    return(r)
  }
  
  
  #Seasonal component (s_t)
  s_t<- function(n){
    
    n <- n
    
    s <- numeric(n)
    
    kappa_s<- rgamma(1,1,0.0001)
    sigma_s <- sqrt(1/kappa_s)
    
    initial_WN <- rnorm(n, mean = 0, sd = sigma_s)
    
    # Simulating the seasonal component
    for (t in 1:n) {
      if (t <= 11) {
        s[t] <- initial_WN[t]
      } else {
        s[t] <- rnorm(1, mean = -sum(s[(t-1):(t-11)]), sd = sigma_s)
      }
    }
    
    return(s)
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
  
  
  # Using the trend, seasonal, spatial components and HMM to simulate logLambda_it
  for (i in 1:ndept){
    logLambda_it.vec <- r_t(time) + s_t(time) + x_it(time)
    logLambda_it <- c(logLambda_it, logLambda_it.vec)
  }
  
  df <- data.frame(logLambda_it, Month, Department)
  return(df)
}


SimulatedData<- SimulateMeningococcalData(94,156)

#e_it<- expectedCases

#SimulatedData$mean<- e_it * exp(SimulatedData$logLambda_it)

#SimulatedData$count<- rpois(length(SimulatedData$mean), SimulatedData$mean)

library(dplyr)
SimulatedData <- SimulatedData %>%
  group_by(Department) %>% mutate(Year=seq(1985,1997.917,1/12))
plot(SimulatedData$Year,SimulatedData$logLambda_it, type="l")



#Spatial component (u_i)
# grid of spatial points for n=94 departments
u_i<- function(n){
  n <- n  
  x <- runif(n)
  y <- runif(n)
  spatial_points <- data.frame(x = x, y = y)
  
  kappa_u <- 0.5
  
  # spatial weights matrix
  W <- matrix(0, n, n)
  
  # computing spatial weights
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        
        dist_ij <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
        
        # spatial weight using inverse Euclidean distance
        W[i, j] <- 1/dist_ij
      }
    }
  }
  
  # Simulating the GIA process
  simulated_component <- matrix(rnorm(n), n, 1)
  for (i in 1:n) {
    sum_wy <- sum(W[i, ] * simulated_component)
    simulated_component[i] <- rnorm(1, mean = sum_wy / sum(W[i, ]), sd = 1 / sqrt((kappa_u*sum(W[i, ]))))
  }
  
  return(simulated_component)
  
}

