packages<- c("ar.matrix","sp", "sf", "mapview", "spdep", "ggplot2", "gganimate","gifski", "transformr", "viridis", "rgdal", "tidyverse", "cmdstanr", "rstan", "bayesplot", "posterior", "mvtnorm", "matrixStats")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))
}))

#Westmidlands.shp <- st_read("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/LAD_MAY_2023_UK_BUC_V2.shp")
#France.shp <- st_read("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/ModifiedFranceShapefile.shp")
#mapview(France.shp)

#ADJACENCY MATRIX
#Westmidlands_adjmat <- nb2mat(poly2nb(Westmidlands.shp), style = "B", zero.policy = TRUE)
#France_adjmat <- nb2mat(poly2nb(France.shp), style = "B", zero.policy = TRUE)

time<- time
ndept<- ndept
#time<- 60
#ndept<- length(diag(Westmidlands_adjmat))
#ndept<- length(diag(France_adjmat))

source("Components.R")

#Markov chain (x_it)
MChain<- function(time){
  
  transition_matrix <- matrix(c(0.8, 0.2, 0.4, 0.6), nrow = 2, byrow = TRUE)
  
  initial_state <- 0  
  
  states <- numeric(time)
  
  states[1] <- initial_state
  for (i in 2:time) {
    current_state <- states[i - 1]
    next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
    states[i] <- next_state
  }
  return(states)
}

SimulateMeningococcalData <- function(Model, adj.matrix, time, e_it=matrix(1000,length(diag(adj.matrix)),time)){
  
  time<- length(r)
  ndept<- length(u)
  
  logLambda_it<- matrix(NA, ndept, time)
  EpidemicIndicator<- matrix(NA, ndept, time)
  
  if(Model == 0){
    EpidemicIndicator<- matrix(0, ndept, time)
    for (i in 1:ndept){
      logLambda_it[i, ] <- r + s
    }
    logLambda_it = logLambda_it + u
    
    y_it= matrix(rpois(ndept*time,e_it*exp(logLambda_it)),ndept,time)
    
    return(list(y_it, EpidemicIndicator)) 
  }
  
  else if(Model == 1 || Model == 2 || Model == 4 || Model == 5){
    B<- 1.09
    
    for (i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time)
      logLambda_it[i, ] <- r + s + EpidemicIndicator[i, ] * B
    }
    
    logLambda_it = logLambda_it + u
    
    y_it= matrix(rpois(ndept*time,e_it*exp(logLambda_it)),ndept,time)
    
    return(list(y_it, EpidemicIndicator)) 
  }
  
  else if(Model == 3 || Model == 6){
    B<- c(0.58, 0.47)
    
    for (i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time)
      logLambda_it[i, ] <- r + s + EpidemicIndicator[i, ] * B[1] + EpidemicIndicator[i, ] * B[2]
    }
    
    logLambda_it = logLambda_it + u
    
    y_it= matrix(rpois(ndept*time,e_it*exp(logLambda_it)),ndept,time)
    
    return(list(y_it, EpidemicIndicator))
  }
}
