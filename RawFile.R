# Space  ==> i=1,...,I.
# Time  ==> t=1,...,T.

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
    r[2] <- r[1]
    
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
time<- 24
ndept<- length(diag(Westmidlands_adjmat))
SimulatedData<- SimulateMeningococcalData(adj.matrix=Westmidlands_adjmat, time=time)

#5th department
plot(SimulatedData[5,],type="l", xlab="Month", ylab="Case")


#Likelihood Computation via forward filtering
Gamma <- matrix(c(0.8, 0.2, 0.4, 0.6), nrow = 2, byrow = TRUE)
init.density <- c(0.5, 0.5)
ndept <- nrow(SimulatedData)  
time <- ncol(SimulatedData)  
nstate <- ncol(Gamma)  
e.it<- 1000

likelihood <- function(y,r,s,u,Gamma,init.density,e.it) {
  forward.probs <- matrix(NA, ndept, nstate)
  
  for (i in 1:ndept) {
    alpha <- init.density * c(dpois(y[i, 1], lambda = e.it*(exp(r[1] + s[1] + u[i]))),
                              dpois(y[i, 1], lambda = e.it*(exp(r[1] + s[1] + u[i] + 1))))
    
    for (t in 2:time) {
      alpha <- alpha %*% Gamma * c(dpois(y[i, t], lambda = e.it*(exp(r[t] + s[t] + u[i]))),
                                   dpois(y[i, t], lambda = e.it*(exp(r[t] + s[t] + u[i] + 1))))
    }
    
    forward.probs[i,] <- alpha
  }
  
  lklhd<- forward.probs[,1] + forward.probs[,2]
  full.likelihood<- prod(lklhd)
  return(list(forward.probs=forward.probs, full.likelihood=full.likelihood))
}

likelihood(y_it,r,s,u,Gamma,init.density,e.it)


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

loglikelihood(y_it,r,s,u,Gamma,init.density,e.it)


##########################################################################################

#logLikelihood Computation via forward filtering - Using LogSumExp function

library(matrixStats)

loglikelihood2 <- function(y,r,s,u,Gamma,init.density,e.it) {
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
  
  #return(list(log.forward.probs=log.forward.probs, full.log.likelihood=sum(full.log.likelihood)))
  return(sum(full.log.likelihood))
}

loglikelihood2(y_it,r,s,u,Gamma,init.density,e.it)


####################################################################################################

#loglikelihood via forward filtering - Scaling method

likelihoodSc <- function(y,r,s,u,Gamma,init.density,e.it) {
  scaledforward.probs <- matrix(NA, ndept, nstate)
  lscaleVec <- numeric(length(ndept))
  for (i in 1:ndept) {
    alpha <- init.density * c(dpois(y[i, 1], lambda = e.it*(exp(r[1] + s[1] + u[i]))),
                              dpois(y[i, 1], lambda = e.it*(exp(r[1] + s[1] + u[i] + 1))))
    lscale <- log(sum(alpha))
    alpha <- alpha/sum(alpha)
    
    
    for (t in 2:time) {
      alpha <- alpha %*% Gamma * c(dpois(y[i, t], lambda = e.it*(exp(r[t] + s[t] + u[i]))),
                                   dpois(y[i, t], lambda = e.it*(exp(r[t] + s[t] + u[i] + 1))))
      lscale <- lscale + log(sum(alpha))
      alpha <- alpha/sum(alpha)
      
    }
    scaledforward.probs[i,] <- alpha
    lscaleVec[i] <- lscale
  }
  
  full.log.likelihood<- sum(lscaleVec)
  return(list(scaledforward.prob=scaledforward.probs, full.log.likelihood=full.log.likelihood))
}

likelihoodSc(y_it,r,s,u,Gamma,init.density,e.it)

#loglikelihood surface contour plot
G12<- seq(0,1, 0.075)

G21<- seq(0,1, 0.075)

gridset<- as.matrix(expand.grid(G12, G21))

G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

z<- matrix(apply(gridset, 1, function(v) loglikelihood2(y_it,r,s,u,G(v[1],v[2]),init.density,e.it)),
           nrow = length(G21)) 
z2<- pmax(z,-5000)
filled.contour(G12,G21,z2, color.palette=terrain.colors, 
               xlab="Prob. of transitioning to hyperendemic period (Gamma_12)",
               ylab="Gamma_21",
               main="Loglikelihood")

#loglikelihood plot, fixing gamma_12
G21<- seq(0,1, 0.005)
MG21<- matrix(G21,length(G21),1)
G<- function(G12,G21){
  G12<- G12
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

z<- matrix(apply(MG21, 1, function(v) loglikelihood2(y_it,r,s,u,G(v[1]),init.density,e.it)$full.log.likelihood),
           nrow = length(G21)) 

df<- data.frame(G21=G21, loglik=z)
plot(df$G21, df$loglik,type="l")

#MCMC for Gamma_12 and Gamma_21

G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

num_iteration<- 3000
MC_chain<- matrix(rep(NA,num_iteration*2), nrow=num_iteration, ncol=2)
MC_chain[1,]<- c(runif(1), runif(1))
acceptedG12<- 0
acceptedG21<- 0

for(i in 2:num_iteration){
  proposedG12<- abs(rnorm(1,mean=MC_chain[i-1,1], sd=0.02))
  proposedG21<- abs(rnorm(1,mean=MC_chain[i-1,2], sd=0.02))
  
  likelihoodcurrent12<- loglikelihood2(y_it,r,s,u,G(MC_chain[i-1,1],MC_chain[i-1,2]),init.density,e.it)
  likelihoodproposed12<- loglikelihood2(y_it,r,s,u,G(proposedG12,MC_chain[i-1,2]),init.density,e.it)
  
  priorcurrent12<- dbeta(MC_chain[i-1,1], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed12<- dbeta(proposedG12, shape1 = 1, shape2 = 1, log=TRUE) 
  
  mh.ratioG12<- exp(likelihoodproposed12 + priorproposed12 - likelihoodcurrent12 - priorproposed12)
  
  if(!is.na(mh.ratioG12) && runif(1) < mh.ratioG12){
    acceptedG12<- acceptedG12 + 1
    MC_chain[i,1]<- proposedG12
  }
  else{
    MC_chain[i,1]<- MC_chain[i-1,1]
    acceptedG12<- acceptedG12 + 0
  }
  
  likelihoodcurrent21<- loglikelihood2(y_it,r,s,u,G(MC_chain[i-1,1],MC_chain[i-1,2]),init.density,e.it) 
  likelihoodproposed21<- loglikelihood2(y_it,r,s,u,G(MC_chain[i-1,1], proposedG21),init.density,e.it)
  
  priorcurrent21<- dbeta(MC_chain[i-1,2], shape1 = 1, shape2 = 1, log=TRUE)
  priorproposed21<- dbeta(proposedG21, shape1 = 1, shape2 = 1, log=TRUE)


  mh.ratioG21<- exp(likelihoodproposed21 + priorproposed21 - likelihoodcurrent21 - priorproposed21)
  
  if(!is.na(mh.ratioG21) && runif(1) < mh.ratioG21){
    acceptedG21<- acceptedG21 + 1
    MC_chain[i,2]<- proposedG21
  }
  else{
    MC_chain[i,2]<- MC_chain[i-1,2]
    acceptedG21<- acceptedG21 + 0
  }
}

plot(c(1:num_iteration),MC_chain[,1],type="l")
plot(c(1:num_iteration),MC_chain[,2],type="l")
tail(MC_chain)
mean(MC_chain[,1])
mean(MC_chain[,2])
acceptedG12/num_iteration
acceptedG21/num_iteration


#####################
#Mapping Simulated data
OverallCasesDF<- data.frame(OverallCases=rowSums(SimulatedData), AreaCode=Westmidlands.shp$LAD23CD, id=c(1:30), AreaName=Westmidlands.shp$LAD23NM)

library(sp)
library(ggplot2)
library(viridis)
library(rgdal)
library(tidyverse)

shapefile<- readOGR("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/LAD_MAY_2023_UK_BUC_V2.shp")
shapefile <- fortify(shapefile, region = 'LAD23CD')

shp_OverallCasesDF<- sp::merge(shapefile,OverallCasesDF,
                        by.x="id",
                        by.y="AreaCode",
                        all.x=F)

shp_OverallCasesDF <- arrange(shp_OverallCasesDF, order)

SimulatedData.map <- ggplot(shp_OverallCasesDF, aes(long, lat, fill = OverallCases, group = group)) +
  geom_polygon(col = "white") +
  scale_fill_viridis(option = "B", direction = -1) +  # Reverse the color scale
  coord_equal() + theme_void() +
  ggtitle('Total simulated cases in 24 months') +
  labs(fill = "Cases")

SimulatedData.map
