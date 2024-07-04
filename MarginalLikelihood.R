source("GeneralLoglikelihood.R")
source("logDensities.R")
library(Rcpp)
sourceCpp("Fastfunctions.cpp")

ModelEvidence<- function(y, e_it, R, Model, z_it, z_it2, posteriorSamples, num_samples = 10000, 
                         RMCMC = F, StanHMC = F){
  #R MCMC samples 
  if(RMCMC){
    if(Model == 0){
      posteriorSamples<- posteriorSamples[,-(c(1,2,(ncol(posteriorSamples)-1): ncol(posteriorSamples)))]
    }else if(Model == 1 || Model == 2 || Model == 4 || Model ==5){
      posteriorSamples<- posteriorSamples[,-ncol(posteriorSamples)] 
    }else if(Model == 3 || Model ==6){
      posteriorSamples<- posteriorSamples
    }
  }else{  
    
    #Stan HMC samples
    posteriorSamples<- posteriorSamples
  }
  
  mu<- colMeans(posteriorSamples)
  varcov<- cov(posteriorSamples)
  
  EvaluatedLogDensities<- numeric(num_samples)
  
  theta <- matrix(nrow = 1, ncol = ncol(posteriorSamples))
  class(theta) <- "numeric" 
  
  for(i in 1:num_samples){
     mvnfast::rmvt(n=1, mu = mu, sigma = varcov, df = 3, A = theta, ncores = 2)
    if(Model==0){
      if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
        EvaluatedLogDensities[i]<- -Inf
      }else{
        EvaluatedLogDensities[i]<- newlogLikelihood(y, e_it, Model, z_it, z_it2, theta) + logPriors(theta, R, Model)
        - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, ncores = 2, log = TRUE) 
      }  
    }else if(Model!=0){
      if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
        EvaluatedLogDensities[i]<- -Inf
      }else{
        EvaluatedLogDensities[i]<- newlogLikelihood(y, e_it, Model, z_it, z_it2, theta) + logPriors(theta, R, Model)
        - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, ncores = 2, log = TRUE) 
      }
    }
  }
  MarginalLikelihood<- log(1) - log(matrixStats::count(EvaluatedLogDensities!= -Inf)) + logSumExp_cpp(EvaluatedLogDensities)
  return(list(MarginalLikelihood, matrixStats::count(EvaluatedLogDensities!= -Inf)))
}