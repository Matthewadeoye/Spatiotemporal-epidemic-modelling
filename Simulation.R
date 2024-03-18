# Space  ==> i=1,...,I.
# Time  ==> t=1,...,T.

packages<- c("ar.matrix","sp", "sf", "mapview", "spdep", "ggplot2", "gganimate","gifski", "transformr", "viridis", "rgdal", "tidyverse", "cmdstanr", "rstan", "bayesplot", "posterior", "mvtnorm", "matrixStats")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))
}))

Westmidlands.shp <- st_read("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/LAD_MAY_2023_UK_BUC_V2.shp")
mapview(Westmidlands.shp)

#ADJACENCY MATRIX
Westmidlands_adjmat <- nb2mat(poly2nb(Westmidlands.shp), style = "B", zero.policy = TRUE)

time<- 60
ndept<- length(diag(Westmidlands_adjmat))

source("Components.R")

SimulateMeningococcalData <- function(adj.matrix, time, e_it=matrix(1000,length(diag(adj.matrix)),time)){
  
  ndept<- length(diag(adj.matrix))
  
  #Markov chain (x_it)
  x_it<- function(num_steps){
    
    transition_matrix <- matrix(c(0.8, 0.2, 0.4, 0.6), nrow = 2, byrow = TRUE)

    num_steps <- num_steps
    initial_state <- 0  
    
    states <- numeric(num_steps)
    
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
  EpidemicIndicator<- matrix(NA, ndept, time)
  for (i in 1:ndept){
    EpidemicIndicator[i,]<- x_it(time)
    logLambda_it[i,] <- r + s + EpidemicIndicator[i,]
  }
  
  #Add spatial component (u_i)
  logLambda_it = logLambda_it + u
  
  y_it= matrix(rpois(ndept*time,e_it*exp(logLambda_it)),ndept,time)
  
  return(list(y_it, EpidemicIndicator))
}

SimulationResults<- SimulateMeningococcalData(adj.matrix=Westmidlands_adjmat, time=time)